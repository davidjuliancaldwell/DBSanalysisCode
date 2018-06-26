function Z=normalize_pp(plv_data,plv_ref)
s_ref1=std(plv_ref')';
Z=(plv_data)./(s_ref1*ones(1,size(plv_data,2)));