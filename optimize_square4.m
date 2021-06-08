clear all;

%%test_ito_branch
%%ensyu2,3

%% Initialization

%結果格納用フォルダ
folder = 'slpi';

if exist(folder, 'dir')
    %error('really remove?');
    rmdir(folder, 's');
end

mkdir(folder);

nside = 12; %一辺の要素数
x0 = 0.3; %許容材料使用量(x0=1.0 -> 100 %), default:0.3
ntry = 100; %最適化試行数

p = 3; %スケーリングのべき乗, default:3
movelim = 0.1; %ムーブリミット
otype = 1; %1:OC, 2:SLP, 3:CONLIN, 4:MMA, 5:NLopt

lx = 100; ly = 100; %各辺の長さ
nx = nside; ny = nside; %各辺の要素数
evol = lx * ly / (nx * ny); %各要素の体積

%% Make mesh

%メッシュの作成(ne,x,y)とスケーリングのない要素剛性行列(ke)
[nnode, nelm, ne, x, y, mprop, free, f, ke] = mesh(lx, ly, nx, ny);

r = x0 * ones(nelm, 1);
x1 = r; x2 = r; l1 = zeros(nelm, 1); u1 = ones(nelm, 1);
object_hist = zeros(ntry, 1);
constr_hist = zeros(ntry, 1);
sens = zeros(nelm, 1);

%% Topology optimization

%最適化ループの開始
for itry = 1:ntry

    %全体剛性行列
    K = zeros(2 * nnode);

    for ie = 1:nelm
        mapn = ne(ie, 1:4);
        map = [2 * mapn - 1 2 * mapn];
        K(map, map) = K(map, map) + ke * r(ie)^p;
    end

    %変位の算出
    u = zeros(2 * nnode, 1);
    u(free) = K(free, free) \ f(free);

    %目的関数・感度の算出
    object(itry) = u(free)' * f(free);
    %%目的関数の比較
    tolerance = 0.0001; % %許容誤差

    if itry >= 2
        O1 = object(itry);
        O0 = object(itry - 1);
        A = abs((O1 - O0) / O1); % %前の目的関数との差

        %%目的関数の差が許容範囲ならループ終了
        if A < tolerance
            break
        end

    end

    for ie = 1:nelm
        mapn = ne(ie, 1:4);
        map = [2 * mapn - 1 2 * mapn];
        ue = u(map);
        sens(ie) = -ue' * (p * r(ie)^(p - 1) * ke) * ue;
    end

    %チェッカーボード避けに上下左右で感度を平均化する
    sens = average_sens(sens, nelm, nx);

    %最適化
    xj = r;

    if (otype == 1)
        r = optimal_criteria(sens, r, nelm, evol, x0, movelim); %OC法
        %     elseif(otype == 2)
        %         r = seq_lin_prog(sens, r, nelm, evol, x0, movelim); %SLP法
        %     elseif(otype == 3)
        %         r = colin_opt(sens, r, nelm, evol, x0, movelim); %CONLIN
        %     elseif(otype == 4)
        %         [r,l1,u1] = mma_opt(sens, xj, nelm, evol, x0, itry, x1, x2, l1, u1, movelim);
        %     elseif(otype == 5)
        %         [r,l1,u1] = lnopt_opt(sens, xj, nelm, evol, x0, itry, x1, x2, l1, u1, movelim);
    end

    x2 = x1; x1 = xj;

    %密度分布のプロット
    figure(1);

    for ie = 1:nelm
        imap = [ne(ie, 1) ne(ie, 2) ne(ie, 3) ne(ie, 4)];
        xl = x(imap); yl = y(imap);
        c = [1 - r(ie) 1 - r(ie) 1 - r(ie)];
        fill(xl, yl, c);
        hold on;
    end

    rfile = sprintf('%s/r%03d.png', folder, itry);
    print(rfile, '-dpng');
    hold off;
    
    Frate=5;
    Frame(itry) = getframe(1);
    v = VideoWriter('topology_animation.avi');
    v.FrameRate = Frate % Framerate
    open(v);
    writeVideo(v,Frame);
    close(v);
    hold off;

    %目的関数のプロット
    figure(2);
    object_hist(itry) = object(itry);
    constr_hist(itry) = sum(r) / nelm;
    plot(1:itry, constr_hist(1:itry)); xlabel('iteration'); ylabel('Volume');

    figure(3);
    plot(1:itry, object_hist(1:itry)); xlabel('iteration'); ylabel('Object function');

end

N = itry % %試行回数の表示

%% Save result

save('-ascii', sprintf('%s/obj.dat', folder), 'object_hist');

figure(2);
plot(1:itry, constr_hist(1:itry)); xlabel('iteration'); ylabel('Volume');
print(sprintf('%s/constraint.png', folder), '-dpng');

figure(3);
plot(1:itry, object_hist(1:itry));
xlabel('iteration'); ylabel('Object function');
print(sprintf('%s/objectfunc.png', folder), '-dpng');

fn = sprintf('%s/%g.txt', folder, object_hist(itry));
disp(fn);
fid = fopen(fn, 'w');
fclose(fid);