w_in = 0.8;             %הגדרת מרחב העבודה
w_out = 10;
h_out = 4;
delta = 0.1;
N = w_out/delta;
M = h_out/delta;
[Y, X] = meshgrid(0:delta:w_out,0:delta:h_out);

phi = zeros(M,N);       %הגדרת תנאי התחלה לאלגוריתם
V1 = 1;
V2 = -1;

phi(:,1) = V1;          %הגדרת תנאי שפה C1
phi(1,:) = V1;
phi(:,end) = V1;
phi(end,:) = V1;

phi(16,46:58)=V2;       %תנאי שפה C2 ת.ז נגמר ב0 וב7
phi(24,46:58)=V2;
phi(16:24,58)=V2;
for i=1:4
    phi(16+i,46-i:58)=V2;
    phi(24-i,46-i:58)=V2;
end

tau=1e-5;               %הפרש מתאים לתנאי התכנסות

for i=1:100000          %אלגוריתם jacobi לחישוב הפוטנציאל
    phi_old = phi;
    for m=2:M-1
        for n=2:N-1
            if phi(m,n)~= V2 %ווידוא אי פגיעה בPEC ובתחום הפנימי
                phi(m,n)=(phi(m-1,n)+phi(m+1,n)+phi(m,n+1)+phi(m,n-1))/4;
            end
        end
    end

    diff_norm = norm(phi - phi_old, 'fro');
    phi_norm = norm(phi, 'fro');
    
    if diff_norm / phi_norm < tau
        fprintf('התכנס לאחר %d איטרציות\n', i);
        break
    end
end

figure;                 %פלטים של פוטנציאל ושדה חשמלי
imagesc(0:delta:w_out, 0:delta:h_out, phi);
axis equal tight
colorbar
xlabel('y [mm]');
ylabel('x [mm]');
title('פוטנציאל V(x,y)');

[Ey, Ex] = gradient(-phi, delta);

figure;
quiver(Y(1:end-1,1:end-1), X(1:end-1,1:end-1), Ey, Ex);
axis equal tight
xlabel('y [mm]');
ylabel('x [mm]');
title('שדה חשמלי E(x,y)');

E_2left=zeros(8);
for i=1:4
    E_2left(i)=sqrt((Ex(16+i,46-i).^2)+(Ey(16+i,46-i).^2));
    E_2left(8-i)=sqrt((Ex(24-i,46-i).^2)+(Ey(24-i,46-i).^2));
end

figure;
hold on;

subplot(2,2,1);
plot((1:100).*delta, (-1)*sqrt((Ex(1,1:100)).^2+(Ey(1,1:100)).^2), 'b');        %מטען משטחתי על C1
xlabel('y [mm]');
ylabel('\eta/\epsilon [q/m^2]');
grid on;
subplot(2,2,2);
plot((1:100).*delta, (-1)*sqrt((Ex(40,1:100)).^2+(Ey(40,1:100)).^2), 'r');
xlabel('y [mm]');
ylabel('\eta/\epsilon [q/m^2]');
grid on;
subplot(2,2,3);
plot((1:40).*delta, (-1)*sqrt((Ex(1:40,1)).^2+(Ey(1:40,1)).^2),   'g');
xlabel('x [mm]');
ylabel('\eta/\epsilon [q/m^2]');
grid on;
subplot(2,2,4);
plot((1:40).*delta, (-1)*sqrt((Ex(1:40,100)).^2+(Ey(1:40,100)).^2),  'm');
xlabel('x [mm]');
ylabel('\eta/\epsilon [q/m^2]');
grid on;

sigma=sum(sqrt((Ex(1,1:100)).^2+(Ey(1,1:100)).^2))+sum(sqrt((Ex(40,1:100)).^2+(Ey(40,1:100)).^2))+sum(sqrt((Ex(1:40,1)).^2+(Ey(1:40,1)).^2))+sum(sqrt((Ex(1:40,100)).^2+(Ey(1:40,100)).^2));
fprintf('צפיפות מטען קווית היא %d', sigma);

figure;
hold on;

subplot(2,2,1);
plot((46:58).*delta, sqrt((Ex(16,46:58)).^2+(Ey(16,46:58)).^2), 'b');        %מטען משטחתי על C2
xlabel('y [mm]');
ylabel('\eta/\epsilon [q/m^2]');
grid on;
subplot(2,2,2);
plot((46:58).*delta, sqrt((Ex(24,46:58)).^2+(Ey(24,46:58)).^2), 'r');
xlabel('y [mm]');
ylabel('\eta/\epsilon [q/m^2]');
grid on;
subplot(2,2,3);
plot((16:24).*delta, sqrt((Ex(16:24,58)).^2+(Ey(16:24,58)).^2),   'g');
xlabel('x [mm]');
ylabel('\eta/\epsilon [q/m^2]');
grid on;
subplot(2,2,4);
plot((1:7).*delta, E_2left(1:7,1),  'm');
xlabel(' [mm]');
ylabel('\eta/\epsilon [q/m^2]');
grid on;