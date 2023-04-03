fn egcd(a: i128, b: i128) -> (i128, i128, i128) {
    if a==0 {
        (b, 0, 1)
    } else {
        let (g, y, x) = egcd(b%a, a);
        (g, x - (b/a) * y, y)
    }
}

pub fn mod_inverse(a : i128, m: i128) -> i128 {
    let (g, x, y) = egcd(a, m);
    assert_eq!(g, 1, "Modular inversion does not exist");
    (x % m)
}

