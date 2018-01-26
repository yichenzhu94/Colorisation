function d = div(px,py) %% compute d = div(p) where p = (px,py) and size(px) == size(py)
% taken from Remy Arbegel tutorial
  [ny,nx] = size(px);
  div_x      =  px - px(:,[1,1:(nx-1)]);
  div_x(:,1) =  px(:,1);
  div_x(:,nx) = -px(:,nx-1);
  div_y      =  py - py([1,1:(ny-1)],:);
  div_y(1,:) =  py(1,:);
  div_y(ny,:) = -py(ny-1,:);
  d = div_x + div_y;
end