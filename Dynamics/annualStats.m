function annualStats(X1, X3, V3rel, X5)

figure;
% color blind safe colors - https://mikemol.github.io/technique/colorblind/2018/02/11/color-safe-palette.html
color_blue = [0	0.45	0.70];
color_orange = [0.90 0.60 0];
color= {color_blue, color_orange};

subplot(4,1,1);
hold on
for k = 1:length(V3rel)   
    h=bar(k,V3rel(k));
    set(h,'FaceColor',color{k});
end
xticklabels([]);
ylabel({'WEC float - spar'; 'Relative velocity'; '(m/s RMS)'});
box on

subplot(4,1,2);
hold on
for k = 1:length(X3)   
    h=bar(k,X3(k));
    set(h,'FaceColor',color{k});
end
xticklabels([]);
ylabel({'FWT platform'; 'Heave (m RMS)'});
box on

subplot(4,1,3);
hold on
for k = 1:length(X1)   
    h=bar(k,X1(k));
    set(h,'FaceColor',color{k});
end
xticklabels([]);
ylabel({'FWT platform'; 'Surge (m RMS)'});
box on

subplot(4,1,4);
hold on
for k = 1:length(X3)   
    h=bar(k,X5(k));
    set(h,'FaceColor',color{k});
end
xticklabels([]);
ylabel({'FWT platform'; 'Pitch (deg RMS)'});
box on

legend('Standalone','Combined FWT-WEC')


% %%
% figure(3);
% subplot(4,1,1);
% % c = categorical({'Baseline FWT','FWT with external heave WEC','FWT with internal heave WEC'});
% P= [0 83 4.2];
% b= bar(P);%,y);
% ylabel('Power (kW)')
% 
% subplot(4,1,3);
% % c = categorical({'Baseline FWT','FWT with external heave WEC','FWT with internal heave WEC'});
% X1= [0.42 0.43 0.40];
% b= bar(X1);%c,y);
% ylabel('Surge (m)')
% 
% subplot(4,1,2);
% % c = categorical({'Baseline FWT','FWT with external heave WEC','FWT with internal heave WEC'});
% X3= [0.08 0.14 0.42];
% b= bar(X3);%c,y);
% ylabel('Heave (m)')
% 
% subplot(4,1,4);
% c = categorical({'Baseline FWT','FWT with external heave WEC','FWT with internal heave WEC'});
% X5= [0.11 0.11 0.06];
% b= bar(c, X5);%c,y);
% ylabel('Pitch (deg)')
% 
% % b(1).CData= 'b';
% % b(2).CData= 'r';
% % b(3).CData= 'g';
