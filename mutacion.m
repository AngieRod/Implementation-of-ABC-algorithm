function [x] = mutacion(x1, x2, u)
      x = x1 + u .* (x1 - x2);
end