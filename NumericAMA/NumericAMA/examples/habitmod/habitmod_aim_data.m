function [param_,np,modname,neq,nlag,nlead,eqname_,eqtype_,endog_,delay_,vtype_] = ...
     habitmod_aim_data()

% habitmod_aim_data()
%     This function will return various information about the AIM model,
%     but will not compute the G and H matrices.

  eqname = cell(12, 1);
  param = cell(9, 1);
  endog = cell(12, 1);
  delay = zeros(12, 1);
  vtype = zeros(12, 1);
  eqtype = zeros(12, 1);

  modname = 'habitmod';
  neq = 12;
  np = 9;
  nlag = 4;
  nlead = 1;

  eqname(1) = cellstr('cnsndt');
  eqname(2) = cellstr('c_');
  eqname(3) = cellstr('yddt');
  eqname(4) = cellstr('rff');
  eqname(5) = cellstr('pdotcpix');
  eqname(6) = cellstr('rreal');
  eqname(7) = cellstr('z');
  eqname(8) = cellstr('P');
  eqname(9) = cellstr('pvdy');
  eqname(10) = cellstr('pvdP');
  eqname(11) = cellstr('pvdz');
  eqname(12) = cellstr('one');
  eqname_ = char(eqname);

  eqtype(1) = 1;     eqtype(2) = 0;     eqtype(3) = 0;   
  eqtype(4) = 0;     eqtype(5) = 0;     eqtype(6) = 1;   
  eqtype(7) = 1;     eqtype(8) = 1;     eqtype(9) = 1;   
  eqtype(10) = 1;     eqtype(11) = 1;     eqtype(12) = 1;   

  eqtype_ = eqtype;

  param(1) = cellstr('c0');
  param(2) = cellstr('lambda');
  param(3) = cellstr('gam');
  param(4) = cellstr('sig');
  param(5) = cellstr('rho');
  param(6) = cellstr('delt');
  param(7) = cellstr('acc');
  param(8) = cellstr('rhoz');
  param(9) = cellstr('bet');
  param_ = char(param);

  endog(1) = cellstr('cnsndt');
  endog(2) = cellstr('c_');
  endog(3) = cellstr('yddt');
  endog(4) = cellstr('rff');
  endog(5) = cellstr('pdotcpix');
  endog(6) = cellstr('rreal');
  endog(7) = cellstr('z');
  endog(8) = cellstr('P');
  endog(9) = cellstr('pvdy');
  endog(10) = cellstr('pvdP');
  endog(11) = cellstr('pvdz');
  endog(12) = cellstr('one');
  endog_ = char(endog);

  delay(1) = 0;     delay(2) = 0;     delay(3) = 0;   
  delay(4) = 0;     delay(5) = 0;     delay(6) = 0;   
  delay(7) = 0;     delay(8) = 0;     delay(9) = 0;   
  delay(10) = 0;     delay(11) = 0;     delay(12) = 0;   

  delay_ = delay;

  vtype(1) = 0;     vtype(2) = 1;     vtype(3) = 0;   
  vtype(4) = 0;     vtype(5) = 0;     vtype(6) = 1;   
  vtype(7) = 1;     vtype(8) = 1;     vtype(9) = 1;   
  vtype(10) = 1;     vtype(11) = 1;     vtype(12) = 2;   

  vtype_ = vtype;



