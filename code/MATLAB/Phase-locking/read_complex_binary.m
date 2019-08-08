function v = read_complex_binary (filename, count)

  % usage: read_complex_binary (filename, [count])
  %
  %  open filename and return the contents as a column vector,
  %  treating them as 32 bit complex numbers
  %

  m = nargchk (1,2,nargin);
  if (m)
    usage (m);
  end

  if (nargin < 2)
    count = Inf;
  end
  
  fprintf('READCOMPLEXBINARY: ');
  
  if exist(filename,'file')==2
      f = fopen(filename, 'rb');
  else
      error('File does not exist, or is not of type file.');
  end
  
  if (f < 0)
    v = 0;
  else
    t = fread (f, [2, count], 'float');
    fclose (f);
    v = t(1,:) + t(2,:)*1i;
    fprintf('\n--> size(v)=%s\n',num2str(size(v)));
    [r, c] = size (v);
    v = reshape (v, c, r);
  end
  fprintf('done\n');
  end