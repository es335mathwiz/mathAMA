% firmvalue_aim_matrices()
%     This script will compute the G and H matrices.

  g = zeros(4, 8);
  h = zeros(4, 12);

  g(17) = g(17) - (1.0*1);
  g(17) = g(17) - (r*1);
  h(33) = h(33) + 1;
  h(37) = h(37) - (-1.0*1);
  g(22) = g(22) + 1;
  g(6) = g(6) - (1.0*1);
  g(6) = g(6) - ((-1.0*delta)*1);
  g(26) = g(26) - 1;
  g(27) = g(27) + 1;
  g(31) = g(31) - (0.0*1);
  g(32) = g(32) + 1;
  g(16) = g(16) - 1;

  cofg = g;
  cofh = h;
