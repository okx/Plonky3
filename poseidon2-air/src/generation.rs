use alloc::vec;
use alloc::vec::Vec;

use num_bigint::BigUint;
use p3_field::PrimeField;
use p3_field::{AbstractField, Field, PrimeField64};
use p3_goldilocks::{Goldilocks, MATRIX_DIAG_16_GOLDILOCKS};
use p3_matrix::dense::RowMajorMatrix;
use p3_maybe_rayon::prelude::*;
use p3_mersenne_31::{from_u62, Mersenne31, POSEIDON2_INTERNAL_MATRIX_DIAG_16_SHIFTS};
use p3_poseidon2::{apply_mat4, M31_RC_16_30_U32};
use core::mem::size_of;
use crate::columns::Poseidon2Cols;
use crate::num_cols;

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

pub fn biguint_to_u64(input: BigUint) -> u64 {
    let digits = input.to_bytes_le();
    let mut byts = [0; 8];
    digits.iter().enumerate().for_each(|(j, x)| byts[j] = *x);
    let x_u64 = u64::from_le_bytes(byts);
    x_u64
}

pub fn biguint_to_u32(input: BigUint) -> u32 {
    let digits = input.to_bytes_le();
    let mut byts = [0; 4];
    digits.iter().enumerate().for_each(|(j, x)| byts[j] = *x);
    let x_u32 = u32::from_le_bytes(byts);
    x_u32
}

#[derive(Debug, PartialEq)]
pub enum FieldType {
    M31,
    GL64,
}

pub fn generate_trace<F: PrimeField, const WIDTH: usize>(
    input: &mut [F; WIDTH],
    round_constants: Vec<[F; WIDTH]>,
    field_type: FieldType,
) -> RowMajorMatrix<F> {
    let n_rows = 32;

    let ncols = num_cols::<WIDTH>();
    let mut trace = RowMajorMatrix::new(vec![F::zero(); n_rows * ncols], ncols);
    let (prefix, rows, suffix) = unsafe { trace.values.align_to_mut::<Poseidon2Cols<F, WIDTH>>() };
    assert!(prefix.is_empty(), "Alignment should match");
    assert!(suffix.is_empty(), "Alignment should match");
    assert_eq!(rows.len(), n_rows);

    rows.iter_mut()
        // .zip(inputs.iter_mut())
        .enumerate()
        .for_each(|(i, row)| {
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
            let is_external_layer = r != 0
                && (((r - 1) < rounds_f_beginning) || (p_end <= (r - 1) && (r - 1) < rounds));

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
                    cols.add_rc[j] = cols.input[j] + round_constants[r - 1][j];
                }
            } else {
                // Mark the selector as internal.
                cols.is_internal = F::one();

                // Apply the round constants only on the first element.
                cols.add_rc.copy_from_slice(&cols.input);
                cols.add_rc[0] = cols.input[0] + round_constants[r - 1][0];
            };

            // Apply the sbox. for all layers
            // TODO: for internla layer, only needs to apply for the first element
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
                    let mut mat4_state = core::array::from_fn(|k| state[j + k]);
                    apply_mat4(&mut mat4_state);
                    state[j..j + 4].clone_from_slice(&mat4_state);
                }
                let sums: [F; 4] = core::array::from_fn(|k| {
                    (0..WIDTH).step_by(4).map(|j| state[j + k]).sum::<F>()
                });
                for j in 0..WIDTH {
                    state[j] += sums[j % 4];
                }
            } else if cols.is_internal == F::one() {
                match field_type {
                    FieldType::M31 => {
                        let part_sum: u64 = state
                            .iter()
                            .skip(1)
                            .map(|x| {
                                let x_u64 = biguint_to_u64(x.as_canonical_biguint());
                                x_u64
                            })
                            .sum();
                        let full_sum = part_sum + (biguint_to_u64(state[0].as_canonical_biguint()));
                        let s0 = part_sum + biguint_to_u64((-state[0]).as_canonical_biguint());
                        state[0] = F::from_canonical_u32(biguint_to_u32(
                            from_u62(s0).as_canonical_biguint(),
                        ));
                        for i in 1..16 {
                            let si = full_sum
                                + ((biguint_to_u64(state[i].as_canonical_biguint()))
                                    << POSEIDON2_INTERNAL_MATRIX_DIAG_16_SHIFTS[i - 1]);
                            state[i] = F::from_canonical_u32(biguint_to_u32(
                                from_u62(si).as_canonical_biguint(),
                            ));
                        }
                    }
                    FieldType::GL64 => {
                        let sum: F = state.iter().cloned().sum();
                        for i in 0..WIDTH {
                            state[i] *= F::from_canonical_u64(Goldilocks::from_f(MATRIX_DIAG_16_GOLDILOCKS[i]).as_canonical_u64());
                            state[i] += sum.clone();
                        }
                    }
                }
            }

            // Copy the state to the output.
            // println!("air, round {}, state {:?}",r, state);
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
