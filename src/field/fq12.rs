use std::ops;
use crate::field::fq12;
use crate::field::fq2::Fq2;
use crate::field::fq6::Fq6;
use crate::field::fq::Fq;

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
            [val.val[2].clone()*self.non_residue.clone(), val.val[0].clone(), val.val[1].clone()],
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
            [self.val[0].clone().forbenius_map(power), self.val[1].clone().forbenius_map(power).mult_by_fq2(self.forbenius_coeffs_c1[power % 12].clone())],
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

#[cfg(test)]
mod tests {
    use std::str::FromStr;
    use num_bigint::BigInt;
    use super::*;

    fn get_fq2(val1: BigInt, val2: BigInt) -> Fq2 {
        let fq_field_modulus = BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208583").unwrap();
        let fq2_non_residue = Fq::fq(fq_field_modulus.clone(), BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap());
        let fq2_forbenius_coeffs_c1 = [Fq::fq(fq_field_modulus.clone(), BigInt::from(1)),Fq::fq(fq_field_modulus.clone(), BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap())];
        return Fq2::fq2(
            [
                Fq::fq(fq_field_modulus.clone(), BigInt::from(val1)),
                Fq::fq(fq_field_modulus.clone(), BigInt::from(val2))
            ],
            fq2_non_residue.clone(),
            fq2_forbenius_coeffs_c1.clone()
        );
    }

    fn get_fq2_str(val1: &str, val2: &str) -> Fq2 {
        let fq_field_modulus = BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208583").unwrap();
        let fq2_non_residue = Fq::fq(fq_field_modulus.clone(), BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap());
        let fq2_forbenius_coeffs_c1 = [Fq::fq(fq_field_modulus.clone(), BigInt::from(1)),Fq::fq(fq_field_modulus.clone(), BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap())];
        return Fq2::fq2(
            [
                Fq::fq(fq_field_modulus.clone(), BigInt::from_str(val1).unwrap()),
                Fq::fq(fq_field_modulus.clone(), BigInt::from_str(val2).unwrap())
            ],
            fq2_non_residue.clone(),
            fq2_forbenius_coeffs_c1.clone()
        );
    }

    fn get_fq6_forbenius_coeff_c1() -> [Fq2; 6] {
        return [
            get_fq2(BigInt::from(1), BigInt::from(0)),
            get_fq2(BigInt::from_str("21575463638280843010398324269430826099269044274347216827212613867836435027261").unwrap(),BigInt::from_str("10307601595873709700152284273816112264069230130616436755625194854815875713954").unwrap()),
            get_fq2(BigInt::from_str("21888242871839275220042445260109153167277707414472061641714758635765020556616").unwrap(), BigInt::from(0)),
            get_fq2(BigInt::from_str("3772000881919853776433695186713858239009073593817195771773381919316419345261").unwrap(), BigInt::from_str("2236595495967245188281701248203181795121068902605861227855261137820944008926").unwrap()),
            get_fq2(BigInt::from_str("2203960485148121921418603742825762020974279258880205651966").unwrap(), BigInt::from(0)),
            get_fq2(BigInt::from_str("18429021223477853657660792034369865839114504446431234726392080002137598044644").unwrap(), BigInt::from_str("9344045779998320333812420223237981029506012124075525679208581902008406485703").unwrap())
        ]
    }

    fn get_fq6_forbenius_coeff_c2() -> [Fq2; 6] {
        return [
            get_fq2(BigInt::from(1), BigInt::from(0)),
            get_fq2(BigInt::from_str("2581911344467009335267311115468803099551665605076196740867805258568234346338").unwrap(),BigInt::from_str("19937756971775647987995932169929341994314640652964949448313374472400716661030").unwrap()),
            get_fq2(BigInt::from_str("2203960485148121921418603742825762020974279258880205651966").unwrap(), BigInt::from(0)),
            get_fq2(BigInt::from_str("5324479202449903542726783395506214481928257762400643279780343368557297135718").unwrap(), BigInt::from_str("16208900380737693084919495127334387981393726419856888799917914180988844123039").unwrap()),
            get_fq2(BigInt::from_str("21888242871839275220042445260109153167277707414472061641714758635765020556616").unwrap(), BigInt::from(0)),
            get_fq2(BigInt::from_str("13981852324922362344252311234282257507216387789820983642040889267519694726527").unwrap(), BigInt::from_str("7629828391165209371577384193250820201684255241773809077146787135900891633097").unwrap())
        ]
    }

    fn get_fq6(val: [&str; 6]) -> Fq6 {
        let fq6_non_residue = get_fq2(BigInt::from(9), BigInt::from(1));
        assert_eq!(val.len(), 6);

        let a = get_fq2(BigInt::from_str(val[0].clone()).unwrap(), BigInt::from_str(val[1].clone()).unwrap());
        let b = get_fq2(BigInt::from_str(val[2].clone()).unwrap(), BigInt::from_str(val[3].clone()).unwrap());
        let c = get_fq2(BigInt::from_str(val[4].clone()).unwrap(), BigInt::from_str(val[5].clone()).unwrap());
        return Fq6::fq6(
            [a.clone(), b.clone(), c.clone()],
            fq6_non_residue,
            get_fq6_forbenius_coeff_c1(),
            get_fq6_forbenius_coeff_c2()
        )
    }

    fn get_fq12_forbenius_coeffs_c1() -> [Fq2; 12] {
        return [
            get_fq2_str("1", "0"),
            get_fq2_str("8376118865763821496583973867626364092589906065868298776909617916018768340080", "16469823323077808223889137241176536799009286646108169935659301613961712198316"),
            get_fq2_str("21888242871839275220042445260109153167277707414472061641714758635765020556617", "0"),
            get_fq2_str("11697423496358154304825782922584725312912383441159505038794027105778954184319", "303847389135065887422783454877609941456349188919719272345083954437860409601"),
            get_fq2_str("21888242871839275220042445260109153167277707414472061641714758635765020556616", "0"),
            get_fq2_str("3321304630594332808241809054958361220322477375291206261884409189760185844239", "5722266937896532885780051958958348231143373700109372999374820235121374419868"),
            get_fq2_str("21888242871839275222246405745257275088696311157297823662689037894645226208582", "0"),
            get_fq2_str("13512124006075453725662431877630910996106405091429524885779419978626457868503", "5418419548761466998357268504080738289687024511189653727029736280683514010267"),
            get_fq2_str("2203960485148121921418603742825762020974279258880205651966", "0"),
            get_fq2_str("10190819375481120917420622822672549775783927716138318623895010788866272024264", "21584395482704209334823622290379665147239961968378104390343953940207365798982"),
            get_fq2_str("2203960485148121921418603742825762020974279258880205651967", "0"),
            get_fq2_str("18566938241244942414004596690298913868373833782006617400804628704885040364344", "16165975933942742336466353786298926857552937457188450663314217659523851788715"),
        ]
    }


    #[test]
    fn test_fq12_add() {
        let fq12_non_residue = get_fq2(BigInt::from(9), BigInt::from(1));
        let fq12_forbenius_coeffs_c1 = get_fq12_forbenius_coeffs_c1();
        let a = get_fq6(["1", "2", "1", "2", "1", "2"]);
        let b = get_fq6(["2", "2", "5", "2", "10", "2"]);
        let a1 = get_fq6(["10", "2", "1", "23", "12", "24"]);
        let b1 = get_fq6(["2", "21", "54", "21", "110", "52"]);
        let c = Fq12::fq12([a, b], fq12_non_residue.clone(), fq12_forbenius_coeffs_c1.clone());
        let c1 = Fq12::fq12([a1, b1], fq12_non_residue.clone(), fq12_forbenius_coeffs_c1.clone());
        let res_a = get_fq6(["11", "4", "2", "25", "13", "26"]);
        let res_b = get_fq6(["4", "23", "59", "23", "120", "54"]);
        let res = Fq12::fq12([res_a, res_b], fq12_non_residue.clone(), fq12_forbenius_coeffs_c1.clone());
        let d = c + c1;
        assert_eq!(d, res);
    }

    #[test]
    fn test_fq12_sub() {
        let fq12_non_residue = get_fq2(BigInt::from(9), BigInt::from(1));
        let fq12_forbenius_coeffs_c1 = get_fq12_forbenius_coeffs_c1();
        let a = get_fq6(["1", "2", "1", "2", "1", "2"]);
        let b = get_fq6(["2", "2", "5", "2", "10", "2"]);
        let a1 = get_fq6(["10", "2", "1", "23", "12", "24"]);
        let b1 = get_fq6(["2", "21", "54", "21", "110", "52"]);
        let c = Fq12::fq12([a, b], fq12_non_residue.clone(), fq12_forbenius_coeffs_c1.clone());
        let c1 = Fq12::fq12([a1, b1], fq12_non_residue.clone(), fq12_forbenius_coeffs_c1.clone());
        let res_a = get_fq6(["9", "0", "0", "21", "11", "22"]);
        let res_b = get_fq6(["0", "19", "49", "19", "100", "50"]);
        let res = Fq12::fq12([res_a, res_b], fq12_non_residue.clone(), fq12_forbenius_coeffs_c1.clone());
        let d = c1 - c;
        assert_eq!(d, res);
    }


    #[test]
    fn test_fq12_mul() {
        let fq12_non_residue = get_fq2(BigInt::from(9), BigInt::from(1));
        let fq12_forbenius_coeffs_c1 = get_fq12_forbenius_coeffs_c1();
        let a = get_fq6(["1", "2", "1", "2", "1", "2"]);
        let b = get_fq6(["2", "2", "5", "2", "10", "2"]);
        let a1 = get_fq6(["10", "2", "1", "23", "12", "24"]);
        let b1 = get_fq6(["2", "21", "54", "21", "110", "52"]);
        let c = Fq12::fq12([a, b], fq12_non_residue.clone(), fq12_forbenius_coeffs_c1.clone());
        let c1 = Fq12::fq12([a1, b1], fq12_non_residue.clone(), fq12_forbenius_coeffs_c1.clone());
        let res_a = get_fq6(["1351", "7679", "7249", "8615", "8183", "8010"]);
        let res_b = get_fq6(["21888242871839275222246405745257275088696311157297823662689037894645226207728", "7036", "140", "5134", "9", "655"]);
        let res = Fq12::fq12([res_a, res_b], fq12_non_residue.clone(), fq12_forbenius_coeffs_c1.clone());
        let d = c * c1;
        assert_eq!(d, res);
    }

    #[test]
    fn test_fq12_inverse() {
        let fq12_non_residue = get_fq2(BigInt::from(9), BigInt::from(1));
        let fq12_forbenius_coeffs_c1 = get_fq12_forbenius_coeffs_c1();
        let a = get_fq6(["1", "2", "1", "2", "1", "2"]);
        let b = get_fq6(["2", "2", "5", "2", "10", "2"]);
        let c = Fq12::fq12([a, b], fq12_non_residue.clone(), fq12_forbenius_coeffs_c1.clone());
        let c_inv = c.inverse();
        let fq12_one = Fq12::fq12([get_fq6(["1", "0", "0", "0", "0", "0"]), get_fq6(["0", "0", "0","0", "0", "0"])],
                                  fq12_non_residue.clone(),
                                  fq12_forbenius_coeffs_c1.clone());
        assert_eq!(c * c_inv, fq12_one);
    }

    #[test]
    fn test_fq12_div() {
        let fq12_non_residue = get_fq2(BigInt::from(9), BigInt::from(1));
        let fq12_forbenius_coeffs_c1 = get_fq12_forbenius_coeffs_c1();
        let a = get_fq6(["1", "2", "1", "2", "1", "2"]);
        let b = get_fq6(["2", "2", "5", "2", "10", "2"]);
        let a1 = get_fq6(["10", "2", "1", "23", "12", "24"]);
        let b1 = get_fq6(["2", "21", "54", "21", "110", "52"]);
        let c = Fq12::fq12([a, b], fq12_non_residue.clone(), fq12_forbenius_coeffs_c1.clone());
        let c1 = Fq12::fq12([a1, b1], fq12_non_residue.clone(), fq12_forbenius_coeffs_c1.clone());
        let d = c.clone() * c1.clone();
        let e = d.clone() / c1.clone();
        assert_eq!(e, c);
    }

    #[test]
    fn test_fq12_forbenius_map() {
        todo!()
    }
}

