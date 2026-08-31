#pragma once

// Lays out the three ROOT 3D polylines (shaft, two arrowhead barbs) that
// draw one dipole moment as an arrow in the visual-mode 3D view.
//
// Extracted verbatim from positronium.cpp (continuing the split of engine,
// experiments and ROOT presentation apart -- see the session notes).
// Textually included at the same point inside positronium.cpp's shared
// anonymous namespace, itself already inside an
// #ifndef POSITRONIUM_VALIDATION_EXECUTABLE region positronium.cpp opens
// well before this #include -- so this header needs no guard of its own.
// Depends on that namespace already having in scope: Vec3, unit(), cross().

void setDipoleArrow(TPolyLine3D& shaft, TPolyLine3D& leftHead, TPolyLine3D& rightHead,
                    const Vec3& centre, const Vec3& dipole) {
    const Vec3 direction = unit(dipole);
    const Vec3 tip = centre + direction * 0.31;
    const Vec3 helper = std::abs(direction.z) < 0.8 ? Vec3{0, 0, 1} : Vec3{0, 1, 0};
    const Vec3 side = unit(cross(direction, helper));
    const Vec3 base = tip - direction * 0.085;
    shaft.SetPoint(0, centre.x, centre.y, centre.z);
    shaft.SetPoint(1, tip.x, tip.y, tip.z);
    leftHead.SetPoint(0, tip.x, tip.y, tip.z);
    leftHead.SetPoint(1, base.x + side.x*0.055, base.y + side.y*0.055, base.z + side.z*0.055);
    rightHead.SetPoint(0, tip.x, tip.y, tip.z);
    rightHead.SetPoint(1, base.x - side.x*0.055, base.y - side.y*0.055, base.z - side.z*0.055);
}
