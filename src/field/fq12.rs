use std::ops;
use crate::field::fq12;
use crate::field::fq2::Fq2;
use crate::field::fq6::Fq6;

#[derive(Debug, Clone)]
pub struct Fq12 {
    non_residue: Fq2,
    val: [Fq6; 2],
    forbenius_coeffs_c1: [Fq2; 12]
}

impl Fq12 {
    pub fn fq12(val: [Fq6; 2], non_residue: Fq2, forbenius_coeffs_c1: [Fq2; 12]) -> Self {
        Self {
            val,
            non_residue,
            forbenius_coeffs_c1
        }
    }

    pub fn mul_by_residue(&self, val: Fq6) -> Fq6 {
        return Fq6::fq6(
            [val.val[2].clone()*self.non_residue, val.val[0].clone(), val.val[1].clone()],
            val.non_residue,
            val.forbenius_coeffs_c1,
            val.forbenius_coeffs_c2
        )
    }

    pub fn inverse(&self) -> Fq12 {
        let t0 = self.val[0].clone() * self.val[0].clone();
        let t1 = self.val[1].clone() * self.val[1].clone();
        let t2 = t0.clone() - self.mul_by_residue(t1.clone());
        let t3 = t2.inverse();
        return Fq12::fq12(
            [(self.val[0].clone() * t3.clone()), -(self.val[1].clone() * t3.clone())],
            self.non_residue.clone(),
            self.forbenius_coeffs_c1.clone()
        )
    }

    pub fn forbenius_map(&self, power: usize) -> Fq12 {
        return Fq12::fq12(
            [self.val[0].clone().forbenius_map(power), self.val[1].clone().forbenius_map(power).mult_by_fq2(self.frobenius_coeffs_c1[power % 12])],
            self.non_residue.clone(),
            self.forbenius_coeffs_c1.clone()
        )
    }

    pub fn cyclotomic_square(&self) -> Fq12 {
        let mut z0 = self.val[0].val[0].clone();
        let mut z4 = self.val[0].val[1].clone();
        let mut z3 = self.val[0].val[2].clone();
        let mut z2 = self.val[1].val[0].clone();
        let mut z1 = self.val[1].val[1].clone();
        let mut z5 = self.val[1].val[2].clone();

        let mut tmp = z0.clone() * z1.clone();

        let fq6_non_residue = self.val[0].non_residue.clone();

        let t0 = (z0.clone()+z1.clone()) * (z0.clone() + fq6_non_residue.clone() * z1.clone()) - tmp.clone() - fq6_non_residue.clone() * tmp.clone();
        let t1 = tmp.clone() + tmp.clone();
        // t2 + t3*y = (z2 + z3*y)^2 = b^2
        tmp = z2.clone() * z3.clone();
        let t2 = (z2.clone() + z3.clone()) * (z2.clone() + fq6_non_residue.clone() * z3.clone()) - tmp.clone() - fq6_non_residue.clone() * tmp.clone();
        let t3 = tmp.clone() + tmp.clone();
        // t4 + t5*y = (z4 + z5*y)^2 = c^2
        tmp = z4.clone() * z5.clone();
        let t4 = (z4.clone() + z5.clone()) * (z4.clone() + fq6_non_residue.clone() * z5.clone()) - tmp.clone() - fq6_non_residue.clone() * tmp.clone();
        let t5 = tmp.clone() + tmp.clone();

        // for A
        // z0 = 3*t0-2*z0
        z0 = t0.clone() - z0.clone();
        z0 = z0.clone() + z0.clone();
        z0 = z0.clone() + t0.clone();
        // z1 = 3 * t1 + 2 * z1
        z1 = t1.clone() + z1.clone();
        z1 = z1.clone() + z1.clone();
        z1 = z1.clone() + t1.clone();

        // for B

        // z2 = 3 * (xi * t5) + 2 * z2
        tmp = fq6_non_residue.clone() * t5.clone();
        z2 = tmp.clone()+z2.clone();
        z2 = z2.clone()+z2.clone();
        z2 = z2.clone()+tmp.clone();

        // z3 = 3 * t4 - 2 * z3
        z3 = t4.clone() - z3.clone();
        z3 = z3.clone() + z3.clone();
        z3 = z3.clone() + t4.clone();

        // for C

        // z4 = 3 * t2 - 2 * z4
        z4 = t2.clone() - z4.clone();
        z4 = z4.clone() + z4.clone();
        z4 = z4.clone() + t2.clone();

        // z5 = 3 * t3 + 2 * z5
        z5 = t3.clone() + z5.clone();
        z5 = z5.clone() + z5.clone();
        z5 = z5.clone() + t3.clone();

        let fa = Fq6::fq6([z0, z4, z3], self.val[0].non_residue.clone(), self.val[0].forbenius_coeffs_c1.clone(), self.val[0].forbenius_coeffs_c2.clone());
        let fb = Fq6::fq6([z2, z1, z5], self.val[0].non_residue.clone(), self.val[0].forbenius_coeffs_c1.clone(), self.val[0].forbenius_coeffs_c2.clone());

        return Fq12::fq12(
            [fa, fb],
            self.non_residue.clone(),
            self.forbenius_coeffs_c1.clone()
        )
    }

    pub fn unitary_inverse(&self) -> Fq12 {
        return Fq12::fq12(
            [self.val[0].clone(), -self.val[1].clone()],
            self.non_residue.clone(),
            self.forbenius_coeffs_c1.clone()
        )
    }

    pub fn cyclotomic_exp(&self, power: usize) -> Fq12 {
        todo!()
    }
}

impl ops::Add<Fq12> for Fq12 {
    type Output = Fq12;
    fn add(self, rhs: Fq12) -> Self::Output {
        return Fq12::fq12(
            [
                self.val[0].clone() + rhs.val[0].clone(),
                self.val[1].clone() + rhs.val[1].clone()
            ],
            self.non_residue,
            self.forbenius_coeffs_c1
        )
    }
}

impl ops::Sub<Fq12> for Fq12 {
    type Output = Fq12;
    fn sub(self, rhs: Fq12) -> Self::Output {
        return Fq12::fq12(
            [
                self.val[0].clone() - rhs.val[0].clone(),
                self.val[1].clone() - rhs.val[1].clone()
            ],
            self.non_residue,
            self.forbenius_coeffs_c1
        )
    }
}

impl ops::Neg for Fq12 {
    type Output = Fq12;
    fn neg(self) -> Self::Output {
        return Fq12::fq12(
            [-self.val[0].clone(), -self.val[1].clone()],
            self.non_residue,
            self.forbenius_coeffs_c1
        )
    }
}

impl PartialEq for Fq12 {
    fn eq(&self, other: &Self) -> bool {
        return self.val[0].clone() == other.val[0].clone()
            && self.val[1].clone() == other.val[1].clone()
    }

    fn ne(&self, other: &Self) -> bool {
        return !(self.val[0].clone() == other.val[0].clone()
            && self.val[1].clone() == other.val[1].clone())
    }
}

impl ops::Mul<Fq12> for Fq12 {
    type Output = Fq12;
    fn mul(self, rhs: Fq12) -> Self::Output {
        let a0 = self.val[0].clone();
        let a1 = self.val[1].clone();
        let b0 = rhs.val[0].clone();
        let b1 = rhs.val[1].clone();
        let v0 = a0.clone() * b0.clone();
        let v1 = a1.clone() * b1.clone();
        let val1 = v0.clone() + self.mul_by_residue(v1.clone());
        let val2 = (a0.clone() + a1.clone()) * (b0.clone() + b1.clone()) - (v0.clone() + v1.clone());

        return Fq12::fq12([val1, val2], self.non_residue, self.forbenius_coeffs_c1);
    }
}

impl ops::Div<Fq12> for Fq12 {
    type Output = Fq12;
    fn div(self, rhs: Fq12) -> Self::Output {
        return self*rhs.inverse();
    }
}

