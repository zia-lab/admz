load energies
e=energies-min(energies);
load exp.pm
a=figure;
hold
s=size(exp);
i=1:s(1);
h=line([1,2],[exp(i),exp(i)]);
set(h,'LineWidth',0.1);
r=size(e);
j=1:r(1);
h=line([3,4],[e(j),e(j)]);
set(h,'LineWidth',0.1);
b=get(a,'CurrentAxes');
set(b,'XTickLabel','')
axis([0 5 0 max(exp)+1000])
set(b,'XTick',[0 6])
%set(b,'YTickLabels',[0 5000 10000 15000 20000 25000 30000 35000 40000 45000])
ylabel('Energy (cm-1)')
