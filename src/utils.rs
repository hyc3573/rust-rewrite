use std::cell::Cell;

pub type Real = f32;

#[macro_export]
macro_rules! pv {
    ($v:expr) => {
        println!("{}: {}", stringify!($v), $v);
    };
    ($v:expr, $($args:tt)*) => {
        println!("{}: {}", stringify!($v), $v);
        pv!($($args)+);
    };
}

pub struct Idgen {
    cnt: Cell<usize>,
}

impl Idgen {
    pub fn new() -> Self {
        Idgen { cnt: Cell::new(0) }
    }

    pub fn get(&self) -> usize {
        let temp = self.cnt.get();
        self.cnt.set(temp + 1);
        temp
    }

    pub fn size(&self) -> usize {
        self.cnt.get()
    }
}
