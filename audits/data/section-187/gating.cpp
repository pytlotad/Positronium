// 187e: isPositronium() is an exact mass equality, so a mass-perturbed
// pair silently loses separationFloor() and the collision boundary.
//
// FIXED in audit 188: the mass test is now banded, so re-running this
// prints "yes" on every row.  Kept as the record of the diagnosis.

