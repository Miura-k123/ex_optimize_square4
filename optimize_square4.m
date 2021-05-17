clear all;

%% Initialization

%���ʊi�[�p�t�H���_
folder = 'slpi';
if exist(folder, 'dir')
    %error('really remove?');
    rmdir(folder, 's');
end
mkdir(folder);


nside = 10; %��ӂ̗v�f��
x0=0.3; %���e�ޗ��g�p��(x0=1.0 -> 100%), default:0.3
ntry = 10; %�œK�����s��

p=3; %�X�P�[�����O�ׂ̂���, default:3
movelim = 0.1; %���[�u���~�b�g
otype = 1; %1:OC, 2:SLP, 3:CONLIN, 4:MMA, 5:NLopt

lx=100; ly=100; %�e�ӂ̒���
nx=nside; ny=nside; %�e�ӂ̗v�f��
evol=lx*ly/(nx*ny); %�e�v�f�̑̐�

%% Make mesh

%���b�V���̍쐬(ne,x,y)�ƃX�P�[�����O�̂Ȃ��v�f�����s��(ke)
[nnode,nelm,ne,x,y,mprop,free,f,ke] = mesh(lx, ly, nx, ny);


r=x0 * ones(nelm,1);
x1= r; x2=r; l1 = zeros(nelm,1); u1=ones(nelm,1);
object_hist = zeros(ntry,1);
constr_hist = zeros(ntry,1);
sens = zeros(nelm,1);

%% Topology optimization

%�œK�����[�v�̊J�n
for itry=1:ntry
    
    %�S�̍����s��
    K=zeros(2*nnode);
    for ie=1:nelm
        mapn=ne(ie,1:4);
        map=[2*mapn-1 2*mapn];
        K(map, map)=K(map, map) + ke*r(ie)^p;
    end
    
    %�ψʂ̎Z�o
    u=zeros(2*nnode, 1);
    u(free)=K(free, free)\f(free);
    
    %�ړI�֐��E���x�̎Z�o
    object=u(free)'*f(free);
    for ie=1:nelm
        mapn=ne(ie, 1:4);
        map = [2*mapn-1 2*mapn];
        ue=u(map);
        sens(ie)=-ue'*(p*r(ie)^(p-1)*ke)*ue;
    end
    
    %�`�F�b�J�[�{�[�h�����ɏ㉺���E�Ŋ��x�𕽋ω�����
    sens = average_sens(sens, nelm, nx);
    
    %�œK��
    xj = r;
    if(otype == 1)
        r = optimal_criteria(sens, r, nelm, evol, x0, movelim); %OC�@
%     elseif(otype == 2)
%         r = seq_lin_prog(sens, r, nelm, evol, x0, movelim); %SLP�@
%     elseif(otype == 3)
%         r = colin_opt(sens, r, nelm, evol, x0, movelim); %CONLIN
%     elseif(otype == 4)
%         [r,l1,u1] = mma_opt(sens, xj, nelm, evol, x0, itry, x1, x2, l1, u1, movelim);
%     elseif(otype == 5)
%         [r,l1,u1] = lnopt_opt(sens, xj, nelm, evol, x0, itry, x1, x2, l1, u1, movelim);
    end
    x2 = x1; x1 = xj;
    
    
    %���x���z�̃v���b�g
    figure(1);
    for ie=1:nelm
        imap=[ne(ie,1) ne(ie,2) ne(ie,3) ne(ie,4) ];
        xl=x(imap) ; yl=y(imap) ;
        c=[1-r(ie) 1-r(ie) 1-r(ie)];
        fill(xl,yl,c);
        hold on;
    end
    rfile = sprintf('%s/r%03d.png', folder, itry);
    print(rfile, '-dpng');
    hold off;
    
    %�ړI�֐��̃v���b�g
    figure(2);
    object_hist(itry) = object;
    constr_hist(itry) = sum(r)/nelm;
    plot(1:itry, constr_hist(1:itry));  xlabel('iteration'); ylabel('Volume');
    
    figure(3);
    plot(1:itry, object_hist(1:itry));  xlabel('iteration'); ylabel('Object function');
end

%% Save result

save('-ascii', sprintf('%s/obj.dat', folder), 'object_hist');

figure(2);
plot(1:itry, constr_hist(1:itry));  xlabel('iteration'); ylabel('Volume');
print(sprintf('%s/constraint.png', folder), '-dpng');


figure(3);
plot(1:itry, object_hist(1:itry));
xlabel('iteration'); ylabel('Object function');
print(sprintf('%s/objectfunc.png', folder), '-dpng');

fn=sprintf('%s/%g.txt', folder, object_hist(itry));
disp(fn);
fid=fopen(fn, 'w');
fclose(fid);
