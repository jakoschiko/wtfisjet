/// The dimension count of the infinitesimal part of a [`Jet`].
///
/// Usually the dimension count is given by the number of variables for which you want
/// to calculate the derivatives for. It determines the size and computation cost of the
/// infinitesimal part. Therefore you should try to keep this number as small as possible.
///
/// Keep in mind that the dimension count must be fixed for a given computation. You cannot
/// mix up jets with different dimension counts.
///
/// [`Jet`]: crate::Jet
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Dim(pub usize);
