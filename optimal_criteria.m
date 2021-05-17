function  r = optimal_criteria(sensitivity, r, nelm, evol, r0, movelim)

eta=0.8;
%２分法の上界，下界初期値の設定
lambda_low= 0;
lambda_high= 1e4;
% lambda_low= min(abs(sensitivity/evol))/100;
% lambda_high=max(abs(sensitivity/evol))*100;


ra=r;      % 更新設計変数の初期値として現設計変数を設定
itry_oc=0;
while(itry_oc<=1000)
    if itry_oc==1000     %2分法を1000回繰り返しても収束しない場合は終了
        disp('Warning; Failed in getting lambda');
        break;
    end
    
    current_vol=0; total_vol=0;
    lambda=(lambda_low+lambda_high)*0.5;
%     lambda=sqrt(lambda_low*lambda_high);
    
    for ie=1:nelm    % 設計変数の更新
        ra_mid=(-sensitivity(ie)/(evol*lambda))^eta*r(ie);
        ra_min=max((1-movelim)*r(ie),1e-4);
        ra_max=min((1+movelim)*r(ie),1);
        
        if ra_mid<=ra_min
            ra(ie)=ra_min;
        elseif ra_mid< ra_max
            ra(ie)=ra_mid;
        else
            ra(ie)=ra_max;
        end
        current_vol=current_vol+evol*ra(ie);    % 材料の総量の計算
        total_vol=total_vol+evol;
    end
    h=current_vol-total_vol*r0;
    
    if abs(h/(total_vol*r0))<0.001 % 十分に指定材料総量に近ければ２分法を終了
        break;
    end
    
    if h < 0
        lambda_high=lambda;     %　２分法区間上限値の更新
    else
        lambda_low=lambda;     %　２分法区間下限値の更新
    end
    itry_oc=itry_oc+1;
end


r=ra;   % 更新された設計変数
