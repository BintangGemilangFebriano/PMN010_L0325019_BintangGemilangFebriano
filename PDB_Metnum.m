
function PDB_MetodeNumerik()

  clc;
  close all;
  format long g;

  % 1. PARAMETER
  x_awal = 0;
  x_akhir = 2;
  y_awal = 1;
  n = 15;

  h = (x_akhir - x_awal) / n;

  f = @(x,y) (4*x.^3 - 2*x + 1) ./ (2+y);

  y_eksak_fungsi = @(x) ...
    -2 + sqrt(2*x.^4 - 2*x.^2 + 2*x + 9);

  x = linspace(x_awal,x_akhir,n+1)';

  y_euler = zeros(n+1,1);
  y_heun  = zeros(n+1,1);
  y_rk4   = zeros(n+1,1);

  y_eksak = y_eksak_fungsi(x);

  y_euler(1) = y_awal;
  y_heun(1) = y_awal;
  y_rk4(1) = y_awal;

  % 2. METODE EULER
  for i = 1:n
    y_euler(i+1) = y_euler(i) + ...
      h*f(x(i),y_euler(i));
  endfor

  % 3. METODE HEUN (RUNGE-KUTTA ORDE 2)
  for i = 1:n

    k1 = f(x(i),y_heun(i));

    y_prediksi = y_heun(i) + h*k1;

    k2 = f(x(i+1),y_prediksi);

    y_heun(i+1) = y_heun(i) + ...
      (h/2)*(k1+k2);

  endfor

  % 4. METODE RUNGE-KUTTA ORDE 4
  for i = 1:n

    xi = x(i);
    yi = y_rk4(i);

    k1 = h*f(xi,yi);

    k2 = h*f(xi+h/2,yi+k1/2);

    k3 = h*f(xi+h/2,yi+k2/2);

    k4 = h*f(xi+h,yi+k3);

    y_rk4(i+1) = yi + ...
      (k1+2*k2+2*k3+k4)/6;

  endfor

  % 5. MENAMPILKAN TABEL ITERASI
  printf('\n============================================\n');
  printf(' PENYELESAIAN PERSAMAAN DIFERENSIAL BIASA\n');
  printf('============================================\n');

  printf('Persamaan : dy/dx = (4*x^3-2*x+1)/(2+y)\n');
  printf('Kondisi awal : y(0) = 1\n');
  printf('Target : y(2)\n');
  printf('Jumlah iterasi : %d\n',n);
  printf('Ukuran langkah : %.12f\n\n',h);

  printf('%2s %10s %15s %15s %15s %15s\n',...
    'i','x_i','Euler','Heun','RK4','Eksak');

  printf('%s\n',repmat('-',1,91));

  for i = 1:n+1

    printf('%2d %10.6f %15.9f %15.9f %15.9f %15.9f\n',...
      i-1,x(i),y_euler(i),y_heun(i),...
      y_rk4(i),y_eksak(i));

  endfor

  % 6. MENGHITUNG GALAT
  nilai_eksak = y_eksak(end);

  galat_euler = abs(nilai_eksak-y_euler(end));
  galat_heun  = abs(nilai_eksak-y_heun(end));
  galat_rk4   = abs(nilai_eksak-y_rk4(end));

  galat_relatif_euler = ...
    galat_euler/abs(nilai_eksak)*100;

  galat_relatif_heun = ...
    galat_heun/abs(nilai_eksak)*100;

  galat_relatif_rk4 = ...
    galat_rk4/abs(nilai_eksak)*100;

  % 7. MENAMPILKAN HASIL AKHIR
  printf('\n================ HASIL AKHIR ================\n');

  printf('Nilai eksak y(2) : %.12f\n\n',nilai_eksak);

  printf('%-10s %16s %16s %18s\n',...
    'Metode','y(2)','Galat absolut','Galat relatif (%)');

  printf('%s\n',repmat('-',1,66));

  printf('%-10s %16.12f %16.12g %18.12g\n',...
    'Euler',y_euler(end),galat_euler,...
    galat_relatif_euler);

  printf('%-10s %16.12f %16.12g %18.12g\n',...
    'Heun',y_heun(end),galat_heun,...
    galat_relatif_heun);

  printf('%-10s %16.12f %16.12g %18.12g\n',...
    'RK4',y_rk4(end),galat_rk4,...
    galat_relatif_rk4);

  % 8. MENYIMPAN DATA KE CSV
  fid = fopen('hasil_iterasi_pdb.csv','w');

  if fid == -1

    warning('File CSV tidak dapat dibuat.');

  else

    fprintf(fid,'iterasi,x,Euler,Heun,RK4,Eksak\n');

    for i = 1:n+1

      fprintf(fid,'%d,%.12f,%.12f,%.12f,%.12f,%.12f\n',...
        i-1,x(i),y_euler(i),y_heun(i),...
        y_rk4(i),y_eksak(i));

    endfor

    fclose(fid);

    printf('\nFile CSV berhasil dibuat.\n');

  endif

  % 9. GRAFIK PERBANDINGAN SOLUSI
  x_halus = linspace(x_awal,x_akhir,400)';

  y_halus = y_eksak_fungsi(x_halus);

  figure('Name','Perbandingan Solusi PDB');

  plot(x_halus,y_halus,'k-','linewidth',2);
  hold on;

  plot(x,y_euler,'b-o','linewidth',1.2);
  plot(x,y_heun,'r-s','linewidth',1.2);
  plot(x,y_rk4,'g-^','linewidth',1.2);

  grid on;

  xlabel('Variabel independen (x)');
  ylabel('Variabel dependen (y)');

  title('Perbandingan Solusi Eksak dan Metode Numerik');

  legend('Solusi eksak','Euler','Heun','Runge-Kutta Orde 4',...
    'location','northwest');

  print('grafik_perbandingan_pdb.png','-dpng','-r150');

  % 10. GRAFIK PERBANDINGAN GALAT
  galat = [galat_euler,galat_heun,galat_rk4];

  figure('Name','Perbandingan Galat Akhir');

  bar(galat);

  set(gca,'yscale','log');

  set(gca,'xtick',1:3,...
    'xticklabel',{'Euler','Heun','RK4'});

  grid on;

  ylabel('Galat absolut (skala logaritmik)');

  title('Perbandingan Galat Absolut pada x = 2');

  print('grafik_galat_pdb.png','-dpng','-r150');

  printf('\nGrafik berhasil dibuat.\n');
  printf('Proses perhitungan selesai.\n');

endfunction

