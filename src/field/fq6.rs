use std::ops;
use num_bigint::BigInt;
use crate::field::fq2::Fq2;
use crate::field::fq::Fq;

// degree ==> 6
// inner_len ==> 2
// self.val.len ==> 6//2 = 3

#[derive(Debug, Clone)]
pub struct Fq6 {
    pub non_residue: Fq2,
    pub val: [Fq2; 3],
    pub forbenius_coeffs_c1: [Fq2; 6],
    pub forbenius_coeffs_c2: [Fq2; 6]
}

impl Fq6 {
    pub fn fq6(val: [Fq2; 3], non_residue: Fq2, forbenius_coeffs_c1: [Fq2; 6], forbenius_coeffs_c2: [Fq2; 6])  -> Self {
        Self {
            val, non_residue, forbenius_coeffs_c1, forbenius_coeffs_c2
        }
    }

    pub fn inverse(&self) -> Self {
        let a0 = self.val[0].clone();
        let a1 = self.val[1].clone();
        let a2 = self.val[2].clone();

        let t0 = a0.clone() * a0.clone();
        let t1 = a1.clone() * a1.clone();
        let t2 = a2.clone() * a2.clone();
        let t3 = a0.clone() * a1.clone();
        let t4 = a0.clone() * a2.clone();
        let t5 = a1.clone() * a2.clone();

        let c0 = t0.clone() - (t5.clone() * self.non_residue.clone());
        let c1 = (t2.clone() * self.non_residue.clone()) - t3.clone();
        let c2 = t1.clone() - t4.clone();

        let t6 = (a0.clone() * c0.clone()) + (self.non_residue.clone() * ((a2.clone()*c1.clone())+(a1.clone()*c2.clone())));
        let t7 = t6.inverse();

        return Self {
            val: [
                t7.clone() * c0.clone(),
                t7.clone() * c1.clone(),
                t7.clone() * c2.clone()
            ],
            non_residue: self.non_residue.clone(),
            forbenius_coeffs_c1: self.forbenius_coeffs_c1.clone(),
            forbenius_coeffs_c2: self.forbenius_coeffs_c2.clone()
        }
    }

    pub fn mult_by_fq2(&self, fq2: Fq2) -> Self {
        return Fq6::fq6(
            [
                fq2.clone() * self.val[0].clone(),
                fq2.clone() * self.val[1].clone(),
                fq2.clone() * self.val[2].clone()
            ],
            self.non_residue.clone(),
            self.forbenius_coeffs_c1.clone(),
            self.forbenius_coeffs_c2.clone()
        )
    }

    pub fn forbenius_map(&self, power: usize) -> Self {
        return Fq6::fq6(
            [
                self.val[0].clone().forbenius_map(power),
                self.forbenius_coeffs_c1[power%6].clone() * self.val[1].clone().forbenius_map(power),
                self.forbenius_coeffs_c2[power%6].clone() * self.val[2].clone().forbenius_map(power)
            ],
            self.non_residue.clone(),
            self.forbenius_coeffs_c1.clone(),
            self.forbenius_coeffs_c2.clone()
        )
    }
}

impl ops::Add<Fq6> for Fq6 {
    type Output = Fq6;
    fn add(self, rhs: Fq6) -> Self::Output {
        Fq6::fq6(
            [
                self.val[0].clone() + rhs.val[0].clone(),
                self.val[1].clone() + rhs.val[1].clone(),
                self.val[2].clone() + rhs.val[2].clone()
            ],
            self.non_residue,
            self.forbenius_coeffs_c1,
            self.forbenius_coeffs_c2
        )
    }
}

impl ops::Sub<Fq6> for Fq6 {
    type Output = Fq6;
    fn sub(self, rhs: Fq6) -> Self::Output {
        Fq6::fq6(
            [
                self.val[0].clone() - rhs.val[0].clone(),
                self.val[1].clone() - rhs.val[1].clone(),
                self.val[2].clone() - rhs.val[2].clone()
            ],
            self.non_residue,
            self.forbenius_coeffs_c1,
            self.forbenius_coeffs_c2
        )
    }
}

impl ops::Mul<Fq6> for Fq6 {
    type Output = Fq6;
    fn mul(self, rhs: Fq6) -> Self::Output {
        let a:[Fq2; 3] = [self.val[0].clone(), self.val[1].clone(), self.val[2].clone()];
        let b:[Fq2; 3] = [rhs.val[0].clone(), rhs.val[1].clone(), rhs.val[2].clone()];
        let v0 = a[0].clone()*b[0].clone();
        let v1 = a[1].clone()*b[1].clone();
        let v2 = a[2].clone()*b[2].clone();

        let x = ((a[1].clone()+a[2].clone())*(b[1].clone()+b[2].clone())-(v1.clone()+v2.clone()))*self.non_residue.clone();
        let y = ((a[0].clone()+a[1].clone())*(b[0].clone()+b[1].clone())-(v0.clone()+v1.clone()));
        let z = (v2.clone() * self.non_residue.clone());
        let w = (((a[0].clone()+a[2].clone())*(b[0].clone()+b[2].clone()))-(v0.clone()+v2.clone()))+v1.clone();
        Fq6::fq6(
            [
                v0.clone() + x.clone(),
                y+z,
                w
            ],
            self.non_residue,
            self.forbenius_coeffs_c1,
            self.forbenius_coeffs_c2
        )
    }
}

impl ops::Div<Fq6> for Fq6 {
    type Output = Fq6;
    fn div(self, rhs: Fq6) -> Self::Output {
        return self*rhs.inverse();
    }
}

impl ops::Neg for Fq6 {
    type Output = Fq6;
    fn neg(self) -> Self::Output {
        return Fq6::fq6([
            -self.val[0].clone(),
            -self.val[1].clone(),
            -self.val[2].clone()
        ],
      self.non_residue.clone(),
self.forbenius_coeffs_c1.clone(),
self.forbenius_coeffs_c2.clone()
        )
    }
}

impl PartialEq for Fq6 {
    fn eq(&self, other: &Self) -> bool {
        return self.val[0].clone() == other.val[0].clone()
            && self.val[1].clone() == other.val[1].clone()
            && self.val[2].clone() == other.val[2].clone()
    }

    fn ne(&self, other: &Self) -> bool {
        return !(self.val[0].clone() == other.val[0].clone()
            && self.val[1].clone() == other.val[1].clone()
            && self.val[2].clone() == other.val[2].clone())
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;
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


    #[test]
    fn test_fq6_add() {
        let a = get_fq6( ["1", "2", "1", "2", "1", "2"]);
        let b = get_fq6(["2", "2", "5", "2", "10", "2"]);
        let res = get_fq6(["3", "4", "6", "4", "11", "4"]);

        let c = a+b;

        assert_eq!(c, res);
    }


    #[test]
    fn test_fq6_sub() {
        let a = get_fq6(["5", "2", "90", "32", "143", "242"]);
        let b = get_fq6(["2", "2", "5", "2", "10", "2"]);
        let res = get_fq6(["3", "0", "85", "30", "133", "240"]);
        let c = a-b;
        assert_eq!(c, res);
    }


    #[test]
    fn test_fq6_mul() {
        let a = get_fq6(["5", "2", "90", "32", "143", "242"]);
        let b = get_fq6(["2", "2", "5", "2", "10", "2"]);
        let res = get_fq6(["7613", "19045", "5945", "25564", "234", "1140"]);
        let c = a*b;
        assert_eq!(c, res);
    }

    #[test]
    fn test_fq6_inverse() {
        let a = get_fq6(["1", "2", "3", "4", "5", "6"]);
        let b = a.inverse();
        let c = a*b;
        let one = get_fq6(["1", "0", "0", "0", "0", "0"]);
        assert_eq!(c, one);
    }

    #[test]
    fn test_fq6_div() {
        let a = get_fq6(["5", "2", "90", "32", "143", "242"]);
        let b = get_fq6(["2", "2", "5", "2", "10", "2"]);
        let c = a.clone()*b.clone();
        let d = c/b;
        assert_eq!(d, a);
    }

    #[test]
    fn test_fq6_forbenius_map() {
        let a = get_fq6(["5", "2", "90", "32", "143", "242"]);
        let b = a.forbenius_map(3);
        let res = get_fq6([
            "5",
            "21888242871839275222246405745257275088696311157297823662689037894645226208581",
            "17062763550631731903611703332118107358157227495574352822570047045133877604628",
            "14924837800098920432735649127671072646516912760483774822158167301823863129239",
            "21758686387092310821963422632510968520708377124222638435959261943560585754933",
            "601372476167518358025151469221647543937875145581712535033851642213170947420"
        ]);
        assert_eq!(b, res);
    }
}