function sensitivity = average_sens(sensitivity, nelm, nx)

    original_sens = sensitivity;

    %�����̏㉺���E�Ŋ��x�𕽋ς���
    for ie=1:nelm
        ie_near = [ie];
        if (ie - nx) > 0 %���[�łȂ���Έ��
            ie_near = [ie_near ie-nx];
        end
        if (ie + nx) <= nelm %��[�łȂ���Έ��
            ie_near = [ie_near ie+nx];
        end
        if mod(ie, nx) ~= 0 %�E�[�łȂ���Έ�E
            ie_near = [ie_near ie+1];
        end
        if mod(ie, nx) ~= 1 %���[�łȂ���Έ��
            ie_near = [ie_near ie-1];
        end
        sensitivity(ie) = mean(original_sens(ie_near));
    end    

end