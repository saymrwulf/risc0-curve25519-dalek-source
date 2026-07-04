// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2021 isis lovecruft
// Copyright (c) 2016-2019 Henry de Valence
// See LICENSE for licensing information.
//
// Authors:
// - isis agora lovecruft <isis@patternsinthevoid.net>
// - Henry de Valence <hdevalence@hdevalence.ca>
#![allow(non_snake_case)]

use core::cmp::Ordering;

use crate::backend::serial::curve_models::{CompletedPoint, ProjectiveNielsPoint, ProjectivePoint};
#[cfg(feature = "precomputed-tables")]
use crate::backend::serial::curve_models::AffineNielsPoint;
#[cfg(feature = "precomputed-tables")]
use crate::window::NafLookupTable8;
use crate::constants;
use crate::edwards::EdwardsPoint;
use crate::scalar::Scalar;
use crate::traits::Identity;
use crate::window::NafLookupTable5;

/// Compute \\(aA + bB\\) in variable time, where \\(B\\) is the Ed25519 basepoint.

/// AENEAS-COMPAT helper: highest index with a nonzero NAF digit (or 0).
fn dsm_top_index(a_naf: &[i8; 256], b_naf: &[i8; 256]) -> usize {
    // AENEAS-COMPAT: always process all 256 digits. Leading zero digits
    // double the identity (a no-op), so the result of `mul` is unchanged;
    // only the variable-time skip is dropped. (The conditional-fold scan
    // fails Aeneas' loop fixed point.)
    let _ = a_naf;
    let _ = b_naf;
    255
}

/// AENEAS-COMPAT helper: the double-and-add main loop over indices i..0,
/// with all borrows rooted in parameters (pure refactor).
#[cfg(feature = "precomputed-tables")]
fn dsm_loop(
    i: usize,
    a_naf: &[i8; 256],
    b_naf: &[i8; 256],
    table_a: &NafLookupTable5<ProjectiveNielsPoint>,
    table_b: &NafLookupTable8<AffineNielsPoint>,
) -> ProjectivePoint {
    let mut r = ProjectivePoint::identity();
    let mut j: usize = i + 1;
    while j > 0 {
        let idx = j - 1;
        let t0 = r.double();
        let t1 = dsm_step_p(t0, table_a, a_naf[idx]);
        let t2 = dsm_step_b(t1, table_b, b_naf[idx]);
        r = t2.as_projective();
        j -= 1;
    }
    r
}

#[cfg(not(feature = "precomputed-tables"))]
fn dsm_loop(
    i: usize,
    a_naf: &[i8; 256],
    b_naf: &[i8; 256],
    table_a: &NafLookupTable5<ProjectiveNielsPoint>,
    table_b: &NafLookupTable5<ProjectiveNielsPoint>,
) -> ProjectivePoint {
    let mut r = ProjectivePoint::identity();
    let mut j: usize = i + 1;
    while j > 0 {
        let idx = j - 1;
        let t0 = r.double();
        let t1 = dsm_step_p(t0, table_a, a_naf[idx]);
        let t2 = dsm_step_b(t1, table_b, b_naf[idx]);
        r = t2.as_projective();
        j -= 1;
    }
    r
}

/// AENEAS-COMPAT helper: one NAF-digit conditional add against a
/// ProjectiveNiels table (hoisted out of the main loop so its body is
/// single-assignment; pure refactor, semantics identical).
fn dsm_step_p(
    t: CompletedPoint,
    table: &NafLookupTable5<ProjectiveNielsPoint>,
    d: i8,
) -> CompletedPoint {
    if d > 0 {
        &t.as_extended() + &table.select(d as usize)
    } else if d < 0 {
        &t.as_extended() - &table.select((-d) as usize)
    } else {
        t
    }
}

/// AENEAS-COMPAT helper: as `dsm_step_p`, for the basepoint table.
#[cfg(feature = "precomputed-tables")]
fn dsm_step_b(
    t: CompletedPoint,
    table: &NafLookupTable8<AffineNielsPoint>,
    d: i8,
) -> CompletedPoint {
    if d > 0 {
        &t.as_extended() + &table.select(d as usize)
    } else if d < 0 {
        &t.as_extended() - &table.select((-d) as usize)
    } else {
        t
    }
}

/// AENEAS-COMPAT helper: as `dsm_step_p` (no precomputed tables).
#[cfg(not(feature = "precomputed-tables"))]
fn dsm_step_b(
    t: CompletedPoint,
    table: &NafLookupTable5<ProjectiveNielsPoint>,
    d: i8,
) -> CompletedPoint {
    dsm_step_p(t, table, d)
}

pub fn mul(a: &Scalar, A: &EdwardsPoint, b: &Scalar) -> EdwardsPoint {
    let a_naf = a.non_adjacent_form(5);

    #[cfg(feature = "precomputed-tables")]
    let b_naf = b.non_adjacent_form(8);
    #[cfg(not(feature = "precomputed-tables"))]
    let b_naf = b.non_adjacent_form(5);

    // Find starting index (AENEAS-COMPAT: helper, upward fold)
    let i = dsm_top_index(&a_naf, &b_naf);

    let table_A = NafLookupTable5::<ProjectiveNielsPoint>::from(A);
    #[cfg(feature = "precomputed-tables")]
    let table_B = &constants::AFFINE_ODD_MULTIPLES_OF_BASEPOINT;
    #[cfg(not(feature = "precomputed-tables"))]
    let table_B =
        &NafLookupTable5::<ProjectiveNielsPoint>::from(&constants::ED25519_BASEPOINT_POINT);

    let r = dsm_loop(i, &a_naf, &b_naf, &table_A, table_B);

    r.as_extended()
}
