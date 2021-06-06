function  r = optimal_criteria(sensitivity, r, nelm, evol, r0, movelim)

eta=0.8;
%�Q���@�̏�E�C���E�����l�̐ݒ�
lambda_low= 0;
lambda_high= 1e4;
% lambda_low= min(abs(sensitivity/evol))/100;
% lambda_high=max(abs(sensitivity/evol))*100;


ra=r;      % �X�V�݌v�ϐ��̏����l�Ƃ��Č��݌v�ϐ���ݒ�
itry_oc=0;
while(itry_oc<=1000)
    if itry_oc==1000     %2���@��1000��J��Ԃ��Ă��������Ȃ��ꍇ�͏I��
        disp('Warning; Failed in getting lambda');
        break;
    end
    
    current_vol=0; total_vol=0;
    lambda=(lambda_low+lambda_high)*0.5;
%     lambda=sqrt(lambda_low*lambda_high);
    
    for ie=1:nelm    % �݌v�ϐ��̍X�V
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
        current_vol=current_vol+evol*ra(ie);    % �ޗ��̑��ʂ̌v�Z
        total_vol=total_vol+evol;
    end
    h=current_vol-total_vol*r0;
    
    if abs(h/(total_vol*r0))<0.0001 % �\���Ɏw��ޗ����ʂɋ߂���΂Q���@���I��
        break;
    end
    
    if h < 0
        lambda_high=lambda;     %�@�Q���@��ԏ���l�̍X�V
    else
        lambda_low=lambda;     %�@�Q���@��ԉ����l�̍X�V
    end
    itry_oc=itry_oc+1;
end


r=ra;   % �X�V���ꂽ�݌v�ϐ�
