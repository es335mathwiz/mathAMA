function [param_,np,modname,neq,nlag,nlead,eqname_,eqtype_,endog_,delay_,vtype_] = ...
     firmvalue_aim_data()

% firmvalue_aim_data()
%     This function will return various information about the AIM model,
%     but will not compute the G and H matrices.

  eqname = cell(4, 1);
  param = cell(2, 1);
  endog = cell(4, 1);
  delay = zeros(4, 1);
  vtype = zeros(4, 1);
  eqtype = zeros(4, 1);

  modname = 'firmvalue';
  neq = 4;
  np = 2;
  nlag = 1;
  nlead = 1;

  eqname(1) = cellstr('VALUE');
  eqname(2) = cellstr('DIVIDEND');
  eqname(3) = cellstr('E_');
  eqname(4) = cellstr('ONE');
  eqname_ = char(eqname);

  eqtype(1) = 1;     eqtype(2) = 1;     eqtype(3) = 0;   
  eqtype(4) = 1;     eqtype_ = eqtype;

  param(1) = cellstr('r');
  param(2) = cellstr('delta');
  param_ = char(param);

  endog(1) = cellstr('v');
  endog(2) = cellstr('div');
  endog(3) = cellstr('e_');
  endog(4) = cellstr('one');
  endog_ = char(endog);

  delay(1) = 0;     delay(2) = 0;     delay(3) = 0;   
  delay(4) = 0;     delay_ = delay;

  vtype(1) = 1;     vtype(2) = 0;     vtype(3) = 1;   
  vtype(4) = 2;     vtype_ = vtype;



