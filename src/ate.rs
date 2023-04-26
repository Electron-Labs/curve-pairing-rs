use crate::field::fq::Fq;
use crate::field::fq2::Fq2;

#[derive(Debug, Clone)]
pub struct AteG1PreComp {
    px: Fq,
    py: Fq
}

#[derive(Debug, Clone)]
pub struct AteEllCoeffs {
    ell0: Fq2,
    ellvw: Fq2,
    ellvv: Fq2
}

#[derive(Debug, Clone)]
pub struct AteG2PreComp {
    qx: Fq2,
    qy: Fq2, //[TODO::]coeffs type??
}

impl AteG1PreComp {
    pub fn new(x: Fq, y: Fq) -> Self {
        Self {
            px: x,
            py: y
        }
    }
}

impl PartialEq for AteG1PreComp {
    fn eq(&self, other: &Self) -> bool {
        self.px == other.px && self.py == other.py
    }

    fn ne(&self, other: &Self) -> bool {
        !(self.px == other.px && self.py == other.py)
    }
}

impl AteEllCoeffs {
    pub fn new(ell0: Fq2, ellvw: Fq2, ellvv: Fq2) -> Self {
        Self {
            ell0,
            ellvw,
            ellvv
        }
    }
}

impl PartialEq for AteEllCoeffs {
    fn eq(&self, other: &Self) -> bool {
        self.ell0 == other.ell0 && self.ellvv == other.ellvv && self.ellvw == other.ellvw
    }

    fn ne(&self, other: &Self) -> bool {
        !(self.ell0 == other.ell0 && self.ellvv == other.ellvv && self.ellvw == other.ellvw)
    }
}

impl AteG2PreComp {
    pub fn new(x: Fq2, y: Fq2) -> Self {
        Self {
            qx: x,
            qy: y
        }
    }
}

impl PartialEq for AteG2PreComp {
    fn eq(&self, other: &Self) -> bool {
        self.qx == other.qx && self.qy == other.qx
    }

    fn ne(&self, other: &Self) -> bool {
        !(self.qx == other.qx && self.qy == other.qx)
    }
}
