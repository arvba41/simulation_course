function dirplotter(f,tval,yval) %  plots direction field using aboved defined (f,tval,yval)
[tm,ym]=meshgrid(tval,yval); % creates a grid of points stored in matrices in the (t,y)-plane
  dt = tval(2) - tval(1); % determines t_spacing and y_spacing based on tval
  dy = yval(2) - yval(1); % determines t_spacing and y_spacing based on yval
  fv = vectorize(f); % prevents matrix algebra; if f is a function handle fv is a string
  if isa(f,'function_handle') % if f is function handle
  fv = eval(fv); % turn fv into a function handle from a string 
  end
  yp=feval(fv,tm,ym); % evaluates differential equation fv from tm and ym from meshgrid
  s = 1./max(1/dt,abs(yp)./dy)*0.35;
  h = ishold; % define h as ishold which returns 1 if hold is on, and 0 if it is off
  quiver(tval,yval,s,s.*yp,0,'.r'); hold on; % plots vectors at each point w direction based on righthand side of differential equation
  quiver(tval,yval,-s,-s.*yp,0,'.r');
  if h % h above defined as ishold
  hold on % retain current plot 
  else % refers to all not h 
  hold off % resets axes properties to default
  end
  axis([tval(1)-dt/2,tval(end)+dt/2,yval(1)-dy/2,yval(end)+dy/2]) % sets the scaling for axes on current graph to ends of tval and yval   
end    