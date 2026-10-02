use skani::chain;
use skani::screen;
use skani::seeding;
use skani::regression;
use skani::params;
use skani::types;

use wasm_bindgen::prelude::*;
use std::io::Cursor;
use needletail::parse_fastx_reader;
use web_sys::console;

fn fasta_string_to_sketch(fasta: &str, sketch_params: &params::SketchParams) -> Option<types::Sketch> {
    let mut new_sketch = types::Sketch::new(
        sketch_params.marker_c,
        sketch_params.c,
        sketch_params.k,
        "wasm_input".to_string(),
        sketch_params.use_aa,
    );
    let cursor = Cursor::new(fasta.as_bytes().to_vec());
    let reader = parse_fastx_reader(cursor);
    if reader.is_err() {
        return None;
    }
    let mut reader = reader.unwrap();
    let mut j = 0;
    while let Some(record) = reader.next() {
        if let Ok(record) = record {
            let seq = record.seq();
            if seq.len() >= params::MIN_LENGTH_CONTIG {
                new_sketch.contigs.push(String::from_utf8(record.id().to_vec()).unwrap());
                new_sketch.contig_lengths.push(seq.len() as types::GnPosition);
                new_sketch.total_sequence_length += seq.len();
                seeding::fmh_seeds(&seq, sketch_params, j as u32, &mut new_sketch, true);
                j += 1;
            }
        }
    }
    if j > 0 { Some(new_sketch) } else { None }
}

#[wasm_bindgen]
pub fn compare_genomes(fasta_a: &str, fasta_b: &str, k: u32, c: u32, marker_c: u32) -> String {
    console_error_panic_hook::set_once();

    let sketch_params = params::SketchParams::new(
        marker_c as usize,
        c as usize,
        k as usize,
        false,
        false,
    );
    let command_params = default_command_params();

    console::log_1(&"Sketching A...".into());
    let ref_sketch = fasta_string_to_sketch(fasta_a, &sketch_params);
    console::log_1(&"Sketching B...".into());
    let query_sketch = fasta_string_to_sketch(fasta_b, &sketch_params);

    if ref_sketch.is_none() || query_sketch.is_none() {
        return "Error: invalid FASTA input".to_string();
    }

    let ref_sketch = ref_sketch.unwrap();
    let query_sketch = query_sketch.unwrap();

    console::log_1(&format!("ref markers: {}, query markers: {}", ref_sketch.marker_seeds.len(), query_sketch.marker_seeds.len()).into());
    let passed_screen = screen::check_markers_quickly(&ref_sketch, &query_sketch, 0.0, true);
    console::log_1(&format!("passed_screen: {}", passed_screen).into());

    let model_opt = regression::get_model(sketch_params.c, true);
    let map_params = chain::map_params_from_sketch(&ref_sketch, false, &command_params, &model_opt);
    let ani_res = chain::chain_seeds(&ref_sketch, &query_sketch, map_params);
    console::log_1(&"chain_seeds done...".into());

    if !ani_res.ani.is_finite() {
        return "{\"error\": \"No reliable alignment found between these genomes (too divergent or insufficient overlap).\"}".to_string();
    }

    format!("{{\"ani\": {}, \"af_ref\": {}, \"af_query\": {}}}", ani_res.ani, ani_res.align_fraction_query, ani_res.align_fraction_ref)
}

fn default_command_params() -> params::CommandParams {
    params::CommandParams {
        screen: false,
        screen_val: 0.0,
        mode: params::Mode::Dist,
        out_file_name: "".to_string(),
        ref_files: vec![],
        query_files: vec![],
        refs_are_sketch: false,
        queries_are_sketch: false,
        robust: false,
        median: false,
        sparse: false,
        full_matrix: false,
        diagonal: false,
        max_results: 10000,
        individual_contig_q: false,
        individual_contig_r: false,
        min_aligned_frac: 0.0,
        both_min_aligned_frac: 0.0,
        keep_refs: false,
        est_ci: false,
        learned_ani: false,
        detailed_out: false,
        distance: false,
        rescue_small: false,
        separate_sketches: false,
        short_header: false,
    }
}

#[wasm_bindgen]
pub fn compare_triangle(
    fastas: Vec<String>,
    names: Vec<String>,
    k: u32,
    c: u32,
    marker_c: u32,
) -> String {
    console_error_panic_hook::set_once();

    if fastas.len() < 2 {
        return "Error: need at least 2 genomes for triangle comparison".to_string();
    }

    let sketch_params = params::SketchParams::new(
        marker_c as usize,
        c as usize,
        k as usize,
        false,
        false,
    );
    let command_params = default_command_params();

    console::log_1(&format!("Sketching {} genomes...", fastas.len()).into());
    let sketches: Vec<Option<types::Sketch>> = fastas
        .iter()
        .enumerate()
        .map(|(i, fasta)| {
            let mut s = fasta_string_to_sketch(fasta, &sketch_params)?;
            s.file_name = names.get(i).cloned().unwrap_or_else(|| format!("genome_{}", i));
            Some(s)
        })
        .collect();

    let sketches: Vec<types::Sketch> = sketches.into_iter().flatten().collect();
    let n = sketches.len();

    if n < 2 {
        return "Error: fewer than 2 genomes could be sketched. Check file format and minimum contig length (500bp).".to_string();
    }

    console::log_1(&format!("{} genomes sketched successfully.", n).into());

    let kmer_to_sketch = screen::kmer_to_sketch_from_refs(&sketches);

    let mut results: Vec<serde_json::Value> = vec![];
    let model_opt = regression::get_model(sketch_params.c, true);

    for i in 0..n - 1 {
        let ref_sketch = &sketches[i];

        let screened = screen::screen_refs(
            0.0,
            &kmer_to_sketch,
            ref_sketch,
            &sketch_params,
            &sketches,
            true,
        );

        for j in screened {
            if j <= i {
                continue;
            }
            let query_sketch = &sketches[j];
            let map_params = chain::map_params_from_sketch(
                ref_sketch,
                false,
                &command_params,
                &model_opt,
            );
            let ani_res = chain::chain_seeds(ref_sketch, query_sketch, map_params);
            if ani_res.ani > 0.0 {
                results.push(serde_json::json!({
                    "i": i,
                    "j": j,
                    "name_i": ref_sketch.file_name,
                    "name_j": query_sketch.file_name,
                    "ani": ani_res.ani,
                    "af_query": ani_res.align_fraction_query,
                    "af_ref": ani_res.align_fraction_ref,
                }));
            }
        }
    }

    console::log_1(&format!("Triangle done — {} pairs compared.", results.len()).into());

    let response = serde_json::json!({
        "n": n,
        "names": names,
        "results": results,
    });

    serde_json::to_string(&response).unwrap()
}
