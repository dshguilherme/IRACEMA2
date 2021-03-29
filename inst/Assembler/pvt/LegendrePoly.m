function coefs = LegendrePoly(n)
if n == 0
  coefs = 1;
elseif n == 1
  coefs = [1 0];
else
  P_nm1 = 1;
  P_n = [1 0];
    for i=1:(n-1);
      P_np1 = ((2*i+1)*[P_n,0] - i*[0,0,P_nm1])/(i+1); % recurrence
      [P_nm1,P_n] = deal(P_n,P_np1); % shift
    end
    coefs = P_np1;
  end