use alloc::vec;
use alloc::vec::Vec;

use p3_field::PrimeField;
use p3_matrix::dense::RowMajorMatrix;
use p3_maybe_rayon::prelude::*;
use p3_mersenne_31::Mersenne31;
use p3_poseidon2::{apply_mat4, RC_16_30_U32};
use tracing::instrument;

use crate::air::matmul_internal;
use crate::columns::Poseidon2Cols;
use crate::{num_cols, MATRIX_DIAG_16_M31_U32};

// // TODO: Take generic iterable
// #[instrument(name = "generate Poseidon2 trace", skip_all)]
// pub fn generate_trace_rows<
//     F: PrimeField,
//     const WIDTH: usize,
//     const SBOX_DEGREE: usize,
//     const SBOX_REGISTERS: usize,
//     const HALF_FULL_ROUNDS: usize,
//     const PARTIAL_ROUNDS: usize,
// >(
//     inputs: Vec<[F; WIDTH]>,
// ) -> RowMajorMatrix<F> {
//     let n = inputs.len();
//     assert!(
//         n.is_power_of_two(),
//         "Callers expected to pad inputs to a power of two"
//     );
//     let ncols = num_cols::<WIDTH, SBOX_DEGREE, SBOX_REGISTERS, HALF_FULL_ROUNDS, PARTIAL_ROUNDS>();
//     println!("ncols: {:?}, WIDTH: {:?}, SBOX_DEGREE: {:?}, SBOX_REGISTERS: {:?}, HALF_FULL_ROUNDS: {:?}, PARTIAL_ROUNDS: {:?}", ncols, WIDTH, SBOX_DEGREE, SBOX_REGISTERS, HALF_FULL_ROUNDS, PARTIAL_ROUNDS);
//     let mut trace = RowMajorMatrix::new(vec![F::zero(); n * ncols], ncols);
//     let (prefix, rows, suffix) = unsafe {
//         trace.values.align_to_mut::<Poseidon2Cols<
//             F,
//             WIDTH,
//             SBOX_DEGREE,
//             SBOX_REGISTERS,
//             HALF_FULL_ROUNDS,
//             PARTIAL_ROUNDS,
//         >>()
//     };
//     assert!(prefix.is_empty(), "Alignment should match");
//     assert!(suffix.is_empty(), "Alignment should match");
//     assert_eq!(rows.len(), n);

//     rows.iter_mut().zip(inputs).for_each(|(row, input)| {
//         generate_trace_rows_for_perm(row, input);
//     });

//     trace
// }

// /// `rows` will normally consist of 24 rows, with an exception for the final row.
// fn generate_trace_rows_for_perm<
//     F: PrimeField,
//     const WIDTH: usize,
//     const SBOX_DEGREE: usize,
//     const SBOX_REGISTERS: usize,
//     const HALF_FULL_ROUNDS: usize,
//     const PARTIAL_ROUNDS: usize,
// >(
//     _row: &mut Poseidon2Cols<
//         F,
//         WIDTH,
//         SBOX_DEGREE,
//         SBOX_REGISTERS,
//         HALF_FULL_ROUNDS,
//         PARTIAL_ROUNDS,
//     >,
//     _input: [F; WIDTH],
// ) {
//     // println!("row input: {:?}, export: {:?}, full_rounds.r0,width0: {:?}, full_rounds.r0,width15: {:?}", _row.inputs, _row.export,
//     // _row.beginning_full_rounds[0].sbox[0].0, _row.beginning_full_rounds[0].sbox[15].0);
//     // println!("_input: {:?}", _input);
// }

pub const NUM_POSEIDON2_COLS: usize = size_of::<Poseidon2Cols<Mersenne31, 16>>();

pub fn generate_trace<
    F: PrimeField,
    const WIDTH: usize,
    // const SBOX_DEGREE: usize,
    // const SBOX_REGISTERS: usize,
    // const HALF_FULL_ROUNDS: usize,
    // const PARTIAL_ROUNDS: usize,
>(
    // _: &ExecutionRecord<F>,
    inputs: &mut Vec<[F; WIDTH]>,
    // _: &mut ExecutionRecord<F>,
) -> RowMajorMatrix<F> {
    let n = inputs.len();
    assert!(
        n.is_power_of_two(),
        "Callers expected to pad inputs to a power of two"
    );
    let ncols = num_cols::<WIDTH>();
    // println!("ncols: {:?}, WIDTH: {:?}, SBOX_DEGREE: {:?}, SBOX_REGISTERS: {:?}, HALF_FULL_ROUNDS: {:?}, PARTIAL_ROUNDS: {:?}", ncols, WIDTH, SBOX_DEGREE, SBOX_REGISTERS, HALF_FULL_ROUNDS, PARTIAL_ROUNDS);
    let mut trace = RowMajorMatrix::new(vec![F::zero(); n * ncols], ncols);
    let (prefix, rows, suffix) = unsafe {
        trace.values.align_to_mut::<Poseidon2Cols<
            F,
            WIDTH,
            // SBOX_DEGREE,
            // SBOX_REGISTERS,
            // HALF_FULL_ROUNDS,
            // PARTIAL_ROUNDS,
        >>()
    };
    assert!(prefix.is_empty(), "Alignment should match");
    assert!(suffix.is_empty(), "Alignment should match");
    assert_eq!(rows.len(), n);

    rows.iter_mut().zip(inputs.iter_mut()).enumerate().for_each(|(i, (row, input))| {
        // generate_trace_rows_for_perm(row, input);
        let cols: &mut Poseidon2Cols<F, WIDTH> = row;
        cols.input = *input;

        let r = i % 31;
        let rounds_f = 8;
        let rounds_p = 22;
        let rounds = rounds_f + rounds_p;
        let rounds_f_beginning = rounds_f / 2;
        let p_end = rounds_f_beginning + rounds_p;

        cols.rounds[r] = F::one();
        let is_initial_layer = r == 0;
        let is_external_layer =
            r != 0 && (((r - 1) < rounds_f_beginning) || (p_end <= (r - 1) && (r - 1) < rounds));

        if is_initial_layer {
            // Mark the selector as initial.
            cols.is_initial = F::one();

            // Don't apply the round constants.
            cols.add_rc.copy_from_slice(&cols.input);
        } else if is_external_layer {
            // Mark the selector as external.
            cols.is_external = F::one();

            // Apply the round constants.
            for j in 0..WIDTH {
                cols.add_rc[j] = cols.input[j] + F::from_wrapped_u32(RC_16_30_U32[r - 1][j]);
            }
        } else {
            // Mark the selector as internal.
            cols.is_internal = F::one();

            // Apply the round constants only on the first element.
            cols.add_rc.copy_from_slice(&cols.input);
            cols.add_rc[0] = cols.input[0] + F::from_wrapped_u32(RC_16_30_U32[r - 1][0]);
        };

        // Apply the sbox. for all layers
        for j in 0..WIDTH {
            cols.sbox_deg_3[j] = cols.add_rc[j] * cols.add_rc[j] * cols.add_rc[j];
            // cols.sbox_deg_7[j] = cols.sbox_deg_3[j] * cols.sbox_deg_3[j] * cols.add_rc[j];
        }

        // What state to use for the linear layer.
        let mut state = if is_initial_layer {
            cols.add_rc
        } else if is_external_layer {
            cols.sbox_deg_3
        } else {
            let mut state = cols.add_rc;
            state[0] = cols.sbox_deg_3[0];
            state
        };

        // Apply either the external or internal linear layer.
        if cols.is_initial == F::one() || cols.is_external == F::one() {
            for j in (0..WIDTH).step_by(4) {
                let mut mat4_state = core::array::from_fn(|k| state[j+k]);
                apply_mat4(&mut mat4_state);
                state[j..j + 4].clone_from_slice(&mat4_state);
            }
            let sums: [F; 4] =
                core::array::from_fn(|k| (0..WIDTH).step_by(4).map(|j| state[j + k]).sum::<F>());
            for j in 0..WIDTH {
                state[j] += sums[j % 4];
            }
        } else if cols.is_internal == F::one() {
            let matmul_constants: [F; WIDTH] = MATRIX_DIAG_16_M31_U32
                .iter()
                .map(|x| F::from_wrapped_u32(*x))
                .collect::<Vec<_>>()
                .try_into()
                .unwrap();
            matmul_internal(&mut state, matmul_constants);
        }

        // Copy the state to the output.
        cols.output.copy_from_slice(&state);

        *input = cols.output;
    });

    // let mut input = [F::one(); WIDTH];
    // let rows = (0..128)
    //     .map(|i| {
    //         let mut row = [F::zero(); NUM_POSEIDON2_COLS];

    //         row
    //     })
    //     .collect::<Vec<_>>();

    // Convert the trace to a row major matrix.
    // let mut trace = RowMajorMatrix::new(
    //     rows.into_iter().flatten().collect::<Vec<_>>(),
    //     NUM_POSEIDON2_COLS,
    // );

    // Pad the trace to a power of two.
    // pad_to_power_of_two::<NUM_POSEIDON2_COLS, F>(&mut trace.values);

    trace
}

pub fn pad_to_power_of_two<const N: usize, T: Clone + Default>(values: &mut Vec<T>) {
    debug_assert!(values.len() % N == 0);
    let mut n_real_rows = values.len() / N;
    if n_real_rows < 16 {
        n_real_rows = 16;
    }
    values.resize(n_real_rows.next_power_of_two() * N, T::default());
}
