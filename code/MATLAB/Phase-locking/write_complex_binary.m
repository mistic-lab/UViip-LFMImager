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
  fprintf('WRITECOMPLEXBINARY: ');
  
  %Make sure that the vector is dimensionally correct
  shape=size(data);
  if shape(1) == 1
      data = data.';
      fprintf('\n--> Flipped vector');
  end

  f = fopen (filename, 'wb');
  
  if (f < 0)
    ret = -1;
  else
    I = real(data);
    Q = imag(data);
    fwrite (f, [I Q].', 'float');
	fprintf("\n--> size(I)=%s \n--> size(Q)=%s \n--> size([I Q].')=%s\n",...,
        num2str(size(I)),num2str(size(Q)),num2str(size([I Q].')));
    ret = fclose (f);
  end
  fprintf('done\n');
end