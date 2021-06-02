function [nnode, nelm, ne, x, y, mprop, free, f, ke] = mesh(lx, ly, nx, ny)

young=2e5; pois=0.27;
mprop = [young pois];


%節点座標
dx=lx/nx;
dy=ly/ny;

x1=0:dx:lx;
y1=0:dy:ly;
x=[]; y=[];
for i=1:ny+1 x=[x x1]; end;
for i=1:nx+1 y=[y; y1;]; end;
y=reshape(y, 1, (nx+1)*(ny+1));

%要素コネクティビティ
for i=1:nx
for j=1:ny
    ie=nx*(j-1)+i;
    ne(ie,1)=(nx+1)*(j-1)+i;
    ne(ie,2)=ne(ie,1)+1;
    ne(ie,3)=ne(ie,2)+nx+1;
    ne(ie,4)=ne(ie,3)-1;
end
end
nnode=(nx+1)*(ny+1);
nelm=nx*ny;

%拘束条件
k=0;
for i=1:ny+1
   in=(nx+1)*(i-1)+1;
   k=k+1; fix(k)=2*in-1;
   k=k+1; fix(k)=2*in;
end

%荷重条件
f=zeros(2*nnode,1);
% in=(nx+1);
% f(2*in)=1;
% ensyu1
in=(nx+1)*(ny/2+1);
f(2*in)=1;

%非拘束リスト
free=setdiff([1:2*nnode],fix);


%要素剛性行列の作成
ke=zeros(8);
xl=[0 dx dx 0];
yl=[0 0 dy dy];
p=1/3^0.5; s=[-p p p -p]; t=[-p -p p p];

for i=1:4
    H2=0.25*[-1+t(i) 1-t(i) 1+t(i) -1-t(i); -1+s(i) -1-s(i) 1+s(i) 1-s(i)];
    ZE=[0 0 0 0;0 0 0 0];
    H3=[H2 ZE; ZE H2];
    dx=H2*xl'; dy=H2*yl';
    J=[dx dy];
    ds=det(J);
    G=inv(J);
    sx=G(1,1); tx=G(1,2); sy=G(2,1); ty=G(2,2);
    H1=[sx tx 0 0; 0 0 sy ty; sy ty sx tx];
    b=H1*H3;
    
    young=mprop(1); pois=mprop(2);
    c=zeros(3);
    c1=young/(1-pois*pois);   c2=pois*c1; c3=young/(2*(1+pois));
    c(1,1)=c1; c(2,2)=c1; c(1,2)=c2; c(2,1)=c2; c(3,3)=c3;
    
    ke=ke+ds*b'*c*b;
end



end