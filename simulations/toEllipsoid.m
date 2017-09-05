function E = toEllipsoid(p)

E = [p(1)   0    0  p(3)   0    0;
    0  p(1)   0    0  p(3)   0;
    0    0  p(2)   0    0  p(4);
    p(3)   0    0  p(5)   0    0;
    0  p(3)   0    0  p(5)   0;
    0    0  p(4)   0    0  p(6)];

end