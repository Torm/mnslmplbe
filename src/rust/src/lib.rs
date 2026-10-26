mod be;
mod iter;
mod util;


use console::Term;
use extendr_api::{extendr, extendr_module, Error, Integers, RMatrix, Result, ToVectorValue};
use extendr_api::prelude::{Rfloat, Rint};
use ndarray::Array2;
use crate::be::{be_progress, BlanketSelect};


// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod mnslmplbe;
    fn mnslmplbe;
}


/// Markov network structure learning > Marginal pseudolikelihood (MPL) >
/// Blanket expansion method
///
/// Apply the blanket expansion MPL method to calculate the blanket membership
/// posterior matrix. Can be used to learn a Markov network structure over
/// categorical variables.
///
/// - domain - The domain of the variables, specified by the number of values
///   each categorical variable can take.
/// - sample - A sample of categorical data, which respects the domain.
/// - maxK - Maximum blanket size.
/// - alpha - Value to use for the BDe prior.
/// - parentSelect - Can be `all` (no filtering), `cl` (combination limit
///   filtering) or `cap` (cap selection).
/// - parentSelectFromK - Do not filter blankets until blankets are of size `k`.
/// - combinationLimit - Limit to use for combination limit filtering.
/// - capFactor - Factor to use for cap selection.
///
/// @export
#[extendr]
#[allow(non_snake_case)]
pub fn mnslmplbe(
    domain: Integers,
    sample: RMatrix<Rint>,
    #[default = "1"] maxK: u32,
    #[default = "0.5"] alpha: f64,
    #[default = "\"all\""] parentSelect: &str,
    #[default = "1"] parentSelectFromK: u32,
    #[default = "16000"] combinationLimit: f64,
    #[default = "0.10"] capFactor: f64,
) -> Result<RMatrix<Rfloat>> {
    let domain = u32s_from_integers(&domain);
    let sample = rmatrix_to_array2_int(&sample);
    let terminal = Term::stderr();
    let criteria = match parentSelect {
        "all" => BlanketSelect::All,
        "cl" | "combinationLimit" => BlanketSelect::CombinationLimit(combinationLimit),
        "cap" => BlanketSelect::Cap(capFactor),
        _ => return Err(Error::Other(String::from("Unknown selection criteria.")))
    };

    let blanket_member_posteriors = be_progress(maxK, parentSelectFromK, &domain, &sample, alpha, criteria, Some(&terminal))?;
    Ok(array2_to_rmatrix_float(&blanket_member_posteriors))
}


fn u32s_from_integers(integers: &Integers) -> Box<[u32]> {
    let mut us = Vec::with_capacity(integers.len());
    for i in integers.iter() {
        us.push(i.to_integer() as u32);
    }
    us.into_boxed_slice()
}


fn rmatrix_to_array2_int(matrix: &RMatrix<Rint>) -> Array2<u32> {
    Array2::from_shape_fn((matrix.nrows(), matrix.ncols()), |(i, j)| matrix[[i, j]].to_integer() as u32)
}


fn array2_to_rmatrix_int(array: &Array2<u32>) -> RMatrix<Rint> {
    RMatrix::new_matrix(array.nrows(), array.ncols(), |i, j| Rint::new(array[(i, j)] as i32))
}


fn rmatrix_to_array2_float(matrix: &RMatrix<Rfloat>) -> Array2<f64> {
    Array2::from_shape_fn((matrix.nrows(), matrix.ncols()), |(i, j)| matrix[[i, j]].to_real())
}


fn array2_to_rmatrix_float(array: &Array2<f64>) -> RMatrix<Rfloat> {
    RMatrix::new_matrix(array.nrows(), array.ncols(), |i, j| Rfloat::from(array[(i, j)]))
}
