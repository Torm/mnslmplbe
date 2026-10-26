//! Iterators over variables, neighbours/blankets and configurations.


/// An iterator over all `k`-combinations of `n` variables.
pub struct CombinationIter {
    n: usize,
    at: Box<[u32]>,
    done: bool,
}


impl CombinationIter {

    pub fn new(n: usize, k: usize) -> Self {
        if n == 0 || k == 0 || k > n {
            return Self { n, at: Box::new([]), done: true }
        }
        let mut variables = Vec::with_capacity(k);
        for i in 0..k-1 {
            variables.push(i as u32);
        };
        variables.push(k.wrapping_sub(2) as u32);
        let variables = variables.into_boxed_slice();
        CombinationIter { n, at: variables, done: false }
    }

    pub fn n(&self) -> usize {
        self.n
    }

    pub fn k(&self) -> usize {
        self.at.len()
    }

    pub fn at(&self) -> &[u32] {
        &self.at
    }

    /// Returns false when all combinations have been iterated through.
    pub fn next(&mut self) -> bool {
        if self.done {
            return false;
        }
        let mut pointer = self.k() - 1;
        loop {
            self.at[pointer] = self.at[pointer].wrapping_add(1);
            if self.at[pointer] as usize == self.n - self.k() + 1 + pointer {
                if pointer != 0 { // All digits are k.
                    pointer -= 1;
                    continue;
                } else {
                    self.done = true;
                    return false;
                };
            } else {
                if pointer == self.k() - 1 {
                    return true;
                } else {
                    pointer += 1;
                    self.at[pointer] = self.at[pointer - 1];
                    continue;
                };
            };
        };
    }

}


#[test]
fn test_combination_iter() {
    test_combination_iter_case(5, 3, 10, "[0 1 2][0 1 3][0 1 4][0 2 3][0 2 4][0 3 4][1 2 3][1 2 4][1 3 4][2 3 4]");
    test_combination_iter_case(4, 3, 4, "[0 1 2][0 1 3][0 2 3][1 2 3]");
    test_combination_iter_case(7, 1, 7, "[0][1][2][3][4][5][6]");
    test_combination_iter_case(1, 1, 1, "[0]");
    test_combination_iter_case(7, 7, 1, "[0 1 2 3 4 5 6]");
    test_combination_iter_case(0, 3, 0, "");
    test_combination_iter_case(3, 5, 0, "");
}


#[cfg(test)]
fn test_combination_iter_case(n: usize, k: usize, vs: usize, eq: &str) {
    let mut iter = CombinationIter::new(n, k);
    let mut c = 0;
    let mut s = String::new();
    while iter.next() {
        c += 1;
        let mut st = format!("[");
        for e in iter.at.iter() {
            st = format!("{}{} ", st, e);
        }
        st = format!("{}", st.trim());
        st = format!("{}]", st);
        s = format!("{}{}", s, st);
    }
    assert_eq!(c, vs);
    assert_eq!(s, eq)
}
