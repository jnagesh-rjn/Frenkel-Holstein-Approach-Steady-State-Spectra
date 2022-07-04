[peaks,locs]=findpeaks(S);
[maxm,maxm_idx]=max(peaks);
if J0<0 & size(locs)>1
 %disp('(I00/I01)*L**2');
 Ncoh(index)=(S(locs(maxm_idx))/S(locs(maxm_idx-1)))*(L^2);
end


