function sensitivity = average_sens(sensitivity, nelm, nx)

    original_sens = sensitivity;

    %自分の上下左右で感度を平均する
    for ie=1:nelm
        ie_near = [ie];
        if (ie - nx) > 0 %下端でなければ一個下
            ie_near = [ie_near ie-nx];
        end
        if (ie + nx) <= nelm %上端でなければ一個上
            ie_near = [ie_near ie+nx];
        end
        if mod(ie, nx) ~= 0 %右端でなければ一個右
            ie_near = [ie_near ie+1];
        end
        if mod(ie, nx) ~= 1 %左端でなければ一個左
            ie_near = [ie_near ie-1];
        end
        sensitivity(ie) = mean(original_sens(ie_near));
    end    

end