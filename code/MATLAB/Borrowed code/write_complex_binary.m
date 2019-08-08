function ret = write_complex_binary (data, filename)

  % usage: write_complex_binary (data, filename)
  %
  %  open filename and write the contents of a complex column vector 
  %  32 bit complex number
  %

  m = nargchk (2,2,nargin);
  if (m)
    usage (m);
  end

  f = fopen (filename, 'wb');
  
  if (f < 0)
    ret = -1;
  else
    I = real(data);
    Q = imag(data);
    disp 'Write'
    size(I)
    size(Q)
    size([I Q].')
    fwrite (f, [I Q].', 'float');
    ret = fclose (f);
  end
end