function varargout = untitled(varargin)
%      H = UNTITLED returns the handle to a new UNTITLED or the handle to
%      the existing singleton*.
%
%      UNTITLED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED.M with the given input arguments.
%
%      UNTITLED('Property','Value',...) creates a new UNTITLED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before untitled is made visible.
function untitled_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled (see VARARGIN)

% Choose default command line output for untitled
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in sel_wave.
function sel_wave_Callback(hObject, ~, ~) 
% hObject    handle to sel_wave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global answer ChosenWave A B C
A=0;
B=0;
C=0;
str = get(hObject,'value');
switch str
    case  1
        ChosenWave = 1;
        answer = inputdlg({'Frequency(Hz)','A(V)','ontime(%)'},'parameters',1,{'100','5','50'});         
    case  2
        ChosenWave = 2;
        answer = inputdlg({'Frequency(Hz)','A(V)','ontime(%)'},'parameters',1,{'100','4','0.5'});
    case  3
        ChosenWave = 3;
        answer = inputdlg({'Frequency(Hz)','A(V)','Phi'},'parameters',1,{'100','5','0'});
    case  4
        ChosenWave = 4;
        answer = inputdlg({'Freq1(Hz)','A1(V)','Phi_1','Freq2?Hz?','A2(V)','Phi_2','Freq_3?Hz?','A3(V)','Phi_3'},'parameters',1,{'50','5','0','100','5','0','150','5','0'});
    case  5
        ChosenWave = 5;
        answer = inputdlg({'InitFreqency(Hz)','T(s)','T?Freqency(Hz)','A(V)'},'parameters',1,{'50','0.1','1500','1'});
    case  6
        ChosenWave = 6;
        answer = inputdlg({'P_Noise(dBm)'},'parameters',1,{'2'});
        
end
% Hints: contents = cellstr(get(hObject,'String')) returns sel_wave contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sel_wave


% --- Executes during object creation, after setting all properties.
function sel_wave_CreateFcn(hObject, ~, ~)
% hObject    handle to sel_wave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(~, ~, ~)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, ~, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% global tired
% axes(handles.axes6);
% tired=imread('tired.jpg');
% imshow(tired);

global ChosenWave answer t a f P f0 fs wave C
cla;
                   
axes(handles.axes1);

switch ChosenWave
    case 1
        a=str2double(answer(2));
        f=str2double(answer(1));
        dr=str2double(answer(3)); 
        if((dr<0)||(dr>100))
        errordlg('Ontime Error.','Error');    
        hDialog=findall(0,'tag','Msgbox_Error Dialog');
        btn_ok=findall(hDialog,'style','pushbutton');
        set(btn_ok,'String','Done');  
        else
        N=10*fs/f+1;
        t = (0:N-1)/fs;
        wave=a*square(2*pi*f*t,dr);
        plot(t,wave);
        axis([0 5/f -1.2*a 1.2*a]);
        xlabel('t/s');
        ylabel('x(t)');
        title('Square TimeDomain');
        grid on;
        end
        
    case 2
        a=str2double(answer(2));
        f=str2double(answer(1));
        dr=str2double(answer(3)); 
        if((dr<0)||(dr>1))
        errordlg('Ontime Error','Error');    
        hDialog=findall(0,'tag','Msgbox_Error Dialog');
        btn_ok=findall(hDialog,'style','pushbutton');
        set(btn_ok,'String','Done');
        else
        N=10*fs/f+1;
        t = (0:N-1)/fs;
        wave=a*sawtooth(2*pi*f*t,dr);
        plot(t,wave);
        xlabel('t/s');
        ylabel('x(t)');
        axis([0 5/f -1.2*a 1.2*a]);
        title('Triangle TimeDomain');
        grid on;
        end
    case 3
        a=str2double(answer(2));
        f=str2double(answer(1));
        phi=str2double(answer(3)); 
        N=10*fs/f+1;
        t = (0:N-1)/fs;
        wave=a*sin(2*pi*f*t+phi);
        plot(t,wave);
        axis([0 5/f -1.2*a 1.2*a]);
        xlabel('t/s');
        ylabel('x(t)');
        grid on;
        title('Sine TimeDomain');
    case 4
        a1=str2double(answer(2));
        f1=str2double(answer(1));
        phi1=str2double(answer(3)); 
        a2=str2double(answer(5));
        f2=str2double(answer(4));
        phi2=str2double(answer(6)); 
        a3=str2double(answer(8));
        f3=str2double(answer(7));
        phi3=str2double(answer(9)); 
        f_0=min(f1,f2);
        f=min(f_0,f3);
        a=a1+a2+a3;
        N=10*fs/f+1;
        t = (0:N-1)/fs;
        x1=a1*sin(2*pi*f1*t+phi1);
        x2=a2*sin(2*pi*f2*t+phi2);
        x3=a3*sin(2*pi*f3*t+phi3);
        wave=x1+x2+x3;
        plot(t,wave);
        axis([0 5/f -1.2*a 1.2*a]);
        xlabel('t/s');
        ylabel('x(t)');
        grid on;
        title('Mutiple Sine TimeDomain');
        
    case 5 
        f0 = str2double(answer(1));
        t1 = str2double(answer(2));
        f1= str2double(answer(3));
        a=str2double(answer(4));
        N=10*fs/f+1;
        t = (0:N-1)/fs;
        wave=a*chirp(t,f0,t1,f1);
        plot(t,wave);
        xlabel('t/s');
        ylabel('x(t)');
        axis ([0 4/f0 -1.2*a 1.2*a]); 
        grid on;
        title('Chirp Signal TimeDoamin');
    case 6
        C=1;
        P = str2double(answer(1));
        N=1000;
        wave=wgn(1,N,P);
        t=(0:999)/fs;
        plot(t,wave);
        xlabel('t/s');
        ylabel('x(t)');
        grid on;
        title('Gauss Noise TimeDomain');
        
end

guidata(hObject, handles);

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, ~, ~)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FirstWave
str = get(hObject,'value');
switch str
    case 1
      FirstWave = 1;
    case 2
      FirstWave = 2;
end



% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, ~, ~)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SecondWave answer1
str = get(hObject,'value');
switch str
    case  1%FIR
        SecondWave = 1;
        answer1 = inputdlg({'Minimum Decrease As(dB)','Maximum Decrease Rp(dB)','fp(Hz)','fst(Hz)'},'parameters',1,{'50','2','1000','2000'});  
    case  2
        SecondWave = 2;
        answer1 = inputdlg({'Minimum Decrease As(dB)','Maximum Decrease Rp(dB)','fst2(Hz)','fst1(Hz)','fp2(Hz)','fp1(Hz)'},'parameters',1,{'55','2','2500','500','2000','1000'});
    case  3
        SecondWave = 3;
        answer1 = inputdlg({'Minimum Decrease As(dB)','Maximum Decrease Rp(dB)','fp(Hz)','fst(Hz)'},'parameters',1,{'50','2','1500','1000'}); 
    case  4
        SecondWave = 4;
        answer1 = inputdlg({'Minimum Decrease As(dB)','Maximum Decrease Rp(dB)','fst2(Hz)','fst1(Hz)','fp2(Hz)','fp1(Hz)'},'parameters',1,{'55','2','2000','1000','2500','500'});
    case  5%IIR
        SecondWave = 5;
        answer1 = inputdlg({'Minimum Decrease As(dB)','Maximum Decrease Rp(dB)','fp(Hz)','fst(Hz)'},'parameters',1,{'50','2','1000','2000'});
    case  6
        SecondWave = 6;
        answer1 = inputdlg({'Minimum Decrease As(dB)','Maximum Decrease Rp(dB)','fst2(Hz)','fst1(Hz)','fp2(Hz)','fp1(Hz)'},'parameters',1,{'55','2','2500','500','2000','1000'});  %带通 
    case  7
        SecondWave = 7;
        answer1 = inputdlg({'Minimum Decrease As(dB)','Maximum Decrease Rp(dB)','f_p(Hz)','f_st(Hz)'},'parameters',1,{'50','2','1500','1000'});
    case  8
        SecondWave = 8;
        answer1 = inputdlg({'Minimum Decrease As(dB)','Maximum Decrease Rp(dB)','fst2(Hz)','fst1(Hz)','fp2(Hz)','fp1(Hz)'},'parameters',1,{'55','2','2000','1000','2500','500'});  %带阻 
        
end


% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(~, ~, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FirstWave wave fs B Fs 

axes(handles.axes2);                    

if(FirstWave==1)  %%Fourier Analyse
    if(B==1)
       N=length(wave);
       wave=wave(1:N);
       n=0:(N-1);
       f=n*Fs/N;
       y=fft(wave);
       W=abs(y);
       % plot(f(1:N/2),abs(W(1:N/2))*2/N);%Attention for abs(y)
       stem(f(1:N/2),abs(W(1:N/2))*2/N);%Attention for abs(y)
       xlim([0,0.5*Fs]);
       xlabel('Analog Frequency(Hz)');
       ylabel('A(V)');
       title('Sound Wave Fourier Analyse');
       grid on;
    else
        N=length(wave);
        y=fft(wave,N); 
        mag=abs(y);    %A
        f=(0:N-1)*fs/N; %real freq
        % plot(f(1:round(N/2)),mag(1:round(N/2))*2/N);
        stem(f(1:round(N/2)),mag(1:round(N/2))*2/N)
        xlim([0,0.5*fs]);
        xlabel('Analog Frequency(Hz)');
        ylabel('Magnitude');
        title('Normal Wave Fourier Analyse');
        axis auto;
        grid on; 
    end
      axes(handles.axes2);   
      
    

else%%P
    if(B==1)  %mp3
     window=boxcar(length(wave));
     N=length(wave); 
     wave=wave(1:N);
     [Pxx,f]=periodogram(wave,window,N,Fs); %Directly
     stem(f,10*log10(Pxx))
     % plot(f,10*log10(Pxx));
     xlim([0,0.5*Fs]);
     xlabel('Analog Frequency(Hz)');
     ylabel('P(dB)');
     title('Soundwave P Analyse');
     grid on;
    else   %Normal Signal
     window=boxcar(length(wave));
     N=length(wave); 
     [Pxx,f]=periodogram(wave,window,N,fs); %Directly
     stem(f,10*log10(Pxx));
     %plot(f,10*log10(Pxx));
     xlim([0,0.5*fs]);
     title('Normal Signal P Analyse');
     xlabel('Analog Frequency(Hz)');
     ylabel('P(dB)');
     grid on;
    end
     
     
 end

       % --- Executes on button press in pushbutton3
       
function pushbutton3_Callback(~, ~, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SecondWave wave answer1 X h n B Fs fs
axes(handles.axes3);


switch SecondWave      
    case 1  
      fhp=str2double(answer1(3));
      rp=str2double(answer1(2));
      as=str2double(answer1(1));
      fhst=str2double(answer1(4)); 
    if(fhp>fhst)
        errordlg('Frequency Error','Error');    
        hDialog=findall(0,'tag','Msgbox_Error Dialog');
        btn_ok=findall(hDialog,'style','pushbutton');
        set(btn_ok,'String','Done');
    else
      
      if as<21
      tr=0.9;
       else if as<44
        tr=3.1;
          else if as<53
            tr=3.3;
          else  as<74;
              tr=5.5;
              end
           end
      end
      
      if(B==1)
          wp=2*pi*fhp/Fs;
          wst=2*pi*fhst/Fs;
      else
          wp=2*pi*fhp/fs;
          wst=2*pi*fhst/fs;
      end
     tr_width=wst-wp;
     N=ceil(2*pi*tr/tr_width)+1;
        n=0:(N-1);
        wc=(wp+wst)/2;
        alpha=(N-1)/2;
        hd=(wc/pi)*sinc((wc/pi)*(n-alpha));
        if as<21
            w_fun=boxcar(N);
        else if as<44
                w_fun=hanning(N);           
            else if as<53
                    w_fun=hamming(N);
                else as<74;
                        w_fun=blackman(N);
                end
            end
        end
        h=hd.*w_fun';        

       [Hr,w1]=zerophase(h);
       Hr=20*log10(abs(Hr));
       plot(w1/pi,Hr);
       axis tight;
       xlabel('Digital Freq.\omega/\pi');
       ylabel('20lg|H(e^j^\omega)|(dB)');
       title('LP Filter TimeDomain');
       grid on;    
       X=filter(h,1,wave);
    end
       
    case 2  %%BP
       fhst=str2double(answer1(3));%
       as=str2double(answer1(1));
       rp=str2double(answer1(2));
       flst=str2double(answer1(4)); 
       fhp=str2double(answer1(5));
       flp=str2double(answer1(6));
     if(~((flst<flp)&&(flp<fhp)&&(fhp<fhst)))
        errordlg('Freq Error','Error');    
        hDialog=findall(0,'tag','Msgbox_Error Dialog');
        btn_ok=findall(hDialog,'style','pushbutton');
        set(btn_ok,'String','Done');
     else
       if as<21
         tr=0.9;
          else if as<44
         tr=3.1;
           else if as<53
            tr=3.3;
           else  as<74;
              tr=5.5;
               end
              end
       end
       if(B==1) 
           whp=2*pi*fhp/Fs;
           whst=2*pi*fhst/Fs;
           wlp=2*pi*flp/Fs;
           wlst=2*pi*flst/Fs;
       else
           whp=2*pi*fhp/fs;
           whst=2*pi*fhst/fs;
           wlp=2*pi*flp/fs;
           wlst=2*pi*flst/fs;
       end
       tr_width1=wlp-wlst;
       N1=ceil(6.2*pi/tr_width1)+1;
       tr_width2=whst-whp;
       N2=ceil(6.2*pi/tr_width2)+1;
       N=max(N1,N2);
       wlc=(wlst+wlp)/2;
       whc=(whst+whp)/2;
       alpha= (N-1)/2;
       n=0:1:N-1;
       m=n-alpha+eps;
       hd=[sin(whc*m)-sin(wlc*m)]./(pi*m);        
       if as<21
           w_bman=boxcar(N);
               else if as<44
                        w_bman=hanning(N);
                    else if as<53
                            w_bman=hamming(N);
                        else as<74;
                                w_bman=blackman(N);
                        end
                    end
                end           
       w_bman=w_bman';
       h=hd.*w_bman;            
       [Hr,w1]=zerophase(h);
       Hr=20*log10(abs(Hr));
       plot(w1/pi,Hr);
       axis tight;xlabel('Digital Freq \omega/\pi');
       ylabel('20lg|H(e^j^\omega)|(dB)');
       title('BPF TimeDomain');
       grid on;
       X=filter(h,1,wave);
     end
    case 3
      flp=str2double(answer1(3));
      rp=str2double(answer1(2));
      as=str2double(answer1(1));
      flst=str2double(answer1(4)); 
    if(flp<flst)
        errordlg('Freq Error','Error');    
        hDialog=findall(0,'tag','Msgbox_Error Dialog');
        btn_ok=findall(hDialog,'style','pushbutton');
        set(btn_ok,'String','Done');
    else
       if as<21
          tr=0.9;
         else if as<44
                 tr=3.1;
              else if as<53
                  tr=3.3;
                   else  as<74;
                            tr=5.5;
                   end
             end
       end
       if(B==1)
         wp=2*pi*flp/Fs;
         wst=2*pi*flst/Fs;
       else
         wp=2*pi*flp/fs;
         wst=2*pi*flst/fs;
       end
         tr_width=wp-wst;
         N=ceil(2*pi*tr/tr_width)+1;
         n=0:(N-1);
         wc=(wp+wst)/2;
         alpha=(N-1)/2;
         n=0:1:N-1;
         m=n-alpha+eps;
         hd=[sin(pi*m)-sin(wc*m)]./(pi*m);
        if as<21
            w_fun=boxcar(N);
        else if as<44
                w_fun=hanning(N);              
            else if as<53
                    w_fun=hamming(N);
                else as<74;
                        w_fun=blackman(N);                      
                end
            end
        end
       h=hd.*w_fun';
       axes(handles.axes3);
       [Hr,w1]=zerophase(h);
       Hr=20*log10(abs(Hr));
       plot(w1/pi,Hr);
       axis tight;
       xlabel('Digital Freq\omega/\pi');
       ylabel('20lg|H(e^j^\omega)|(dB)');
       title('HPL TimeDomain');
       grid on;
       X=filter(h,1,wave);
    end
    case 4
       fst2=str2double(answer1(3));
       As=str2double(answer1(1));
       Rp=str2double(answer1(2));
       fst1=str2double(answer1(4)); 
       fp2=str2double(answer1(5));
       fp1=str2double(answer1(6));
     if(~((fp1<fst1)&&(fst1<fst2)&&(fst2<fp2)))
        errordlg('Freq Error','Error');    
        hDialog=findall(0,'tag','Msgbox_Error Dialog');
        btn_ok=findall(hDialog,'style','pushbutton');
        set(btn_ok,'String','Done');
     else  
         if(B==1)
           wp2=2*pi*fp2/Fs;
           wp1=2*pi*fp1/Fs;
           ws2=2*pi*fst2/Fs;
           ws1=2*pi*fst1/Fs;
         else
           wp2=2*pi*fp2/fs;
           wp1=2*pi*fp1/fs;
           ws2=2*pi*fst2/fs;
           ws1=2*pi*fst1/fs;
         end
           As0=50;Rp0=1;
           deltaw=min((ws1-wp1),(wp2-ws2));
           N=ceil((As-7.95)/(2.286*deltaw))+1;
           beta=0.1102*(As-8.7);
           while(As0<As)||(Rp0>Rp)
           n=[0:N-1];
           wd=(kaiser(N,beta))';
           wc1=(ws1+wp1)/2;wc2=(ws2+wp2)/2;
           hd=ideal_lp1(N)-ideal_lp(wc2,N)+ideal_lp(wc1,N);
           h=hd.*wd;
           [db,~,~,~,w]=freqz_m(h,[1]);
           dw=2*pi/1000;
           Rp0=-min(db(round(wp2/dw+1):round(501)));
           As0=-max(db(ws1/dw+1:ws2/dw+1));
           N=N+1;
           end
           plot(w/pi,db);
           xlabel('\omega/\pi');
           ylabel('20lg|H(e^j^\omega)|(dB)');axis([0,1,-90,5]);grid;
           title('BandStop Filter TimeDomain');
           grid on;
           X=filter(h,1,wave);     
     end
    case 5
      fhp=str2double(answer1(3));
      rp=str2double(answer1(2));
      as=str2double(answer1(1));
      fhst=str2double(answer1(4)); 
    if(fhp>fhst)
        errordlg('Freq Error','Error');    
        hDialog=findall(0,'tag','Msgbox_Error Dialog');
        btn_ok=findall(hDialog,'style','pushbutton');
        set(btn_ok,'String','Done');
    else
        if(B==1)
          wp=2*pi*fhp/Fs;
          wst=2*pi*fhst/Fs;
        else
          wp=2*pi*fhp/fs;
          wst=2*pi*fhst/fs;
        end
          wp=wp/pi;
          wst=wst/pi;
          [l,wc] = buttord(wp, wst, rp, as);
          [b,a]=butter(l,wc);
          [H,w]=freqz(b,a);
          H=20*log10(abs(H));
          plot(w/pi,H);
          xlabel('\omega/\pi');ylabel('dB');
          title('LowPass Filter ');
          grid on;
          X=filter(b,a,wave);
    end
    case 6
       fhst=str2double(answer1(3));
       as=str2double(answer1(1));
       rp=str2double(answer1(2));
       flst=str2double(answer1(4)); 
       fhp=str2double(answer1(5));
       flp=str2double(answer1(6));
       if(~((flst<flp)&&(flp<fhp)&&(fhp<fhst)))
        errordlg('Freq Error','Error');    
        hDialog=findall(0,'tag','Msgbox_Error Dialog');
        btn_ok=findall(hDialog,'style','pushbutton');
        set(btn_ok,'String','Done');
       else
           if(B==1)
           whp=2*pi*fhp/Fs;
           whst=2*pi*fhst/Fs;
           wlp=2*pi*flp/Fs;
           wlst=2*pi*flst/Fs;
          else
           whp=2*pi*fhp/fs;
           whst=2*pi*fhst/fs;
           wlp=2*pi*flp/fs;
           wlst=2*pi*flst/fs;
           end
          wlp=wlp/pi;
          wlst=wlst/pi;
          whp=whp/pi;
          whst=whst/pi;
          wp=[wlp,whp];
          wst=[wlst,whst];
          [l,wc] = buttord(wp,wst,rp,as);
          [b,a]=butter(l,wc);
          [H,w]=freqz(b,a);
          H=20*log10(abs(H));
          plot(w/pi,H);
          xlabel('\omega/\pi');ylabel('dB');
          title('BandPass Filter FreqDomain');
          grid on;
          X=filter(b,a,wave);
       end
    case 7
        flp=str2double(answer1(3));
        rp=str2double(answer1(2));
        as=str2double(answer1(1));
        flst=str2double(answer1(4)); 
       if(flp<flst)
        errordlg('Freq Error','Error');    
        hDialog=findall(0,'tag','Msgbox_Error Dialog');
        btn_ok=findall(hDialog,'style','pushbutton');
        set(btn_ok,'String','Done');
       else
           if(B==1)
           wp=2*pi*flp/Fs;
           wst=2*pi*flst/Fs;
           else
           wp=2*pi*flp/fs;
           wst=2*pi*flst/fs;
           end
           wp=wp/pi;
           wst=wst/pi;
           [l,wc] = buttord(wp, wst, rp, as);
           [b,a]=butter(l,wc,'high');
           [H,w]=freqz(b,a);
           H=20*log10(abs(H));
           plot(w/pi,H);
           xlabel('\omega/\pi');ylabel('dB');
           title('HighPass Filter FreqDomain');
           grid on;
           X=filter(b,a,wave);
       end
    case 8
        fhst=str2double(answer1(3));
        as=str2double(answer1(1));
        rp=str2double(answer1(2));
        flst=str2double(answer1(4)); 
        fhp=str2double(answer1(5));
        flp=str2double(answer1(6));
        if(~((flp<flst)&&(flst<fhst)&&(fhst<fhp)))
            errordlg('Freq Error','Error');
            hDialog=findall(0,'tag','Msgbox_Error Dialog');
            btn_ok=findall(hDialog,'style','pushbutton');
            set(btn_ok,'String','Done');
        else  
            if(B==1)
                whp=2*pi*fhp/Fs;
                whst=2*pi*fhst/Fs;
                wlp=2*pi*flp/Fs;
                wlst=2*pi*flst/Fs;
            else
                whp=2*pi*fhp/fs;
                whst=2*pi*fhst/fs;
                wlp=2*pi*flp/fs;
                wlst=2*pi*flst/fs;
           end
          wlp=wlp/pi;
          wlst=wlst/pi;
          whp=whp/pi;
          whst=whst/pi;
          wp=[wlp,whp];
          wst=[wlst,whst];
          [l,wc] = buttord(wp,wst,rp,as);
          [b,a]=butter(l,wc,'stop');
          [H,w]=freqz(b,a);
          H=20*log10(abs(H));
          plot(w/pi,H);
          xlabel('\omega/\pi');ylabel('dB');
          title('BandStop Filter FreqDomain');
          grid on;
          X=filter(b,a,wave);
       end
   
 end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(~, ~, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) ThirdWave
global X Fs B f fs a D
axes(handles.axes4);
D=1;
N=size(X);
N=N(2);
if(B==1)
t=0:(1/Fs):((N-1)/Fs);
plot(t,X);
title('After Filter(TimeDomain)');
grid on;
xlabel('t/s');
ylabel('Magnitude');
sound(X,Fs);
else
    t=0:1/fs:((N-1)/fs);
    plot(t,X);
    title('After Filter(TimeDomain)');
    xlim([0,10/f]);
    ylim([-1.2*a,1.2*a]);
grid on;
xlabel('t/s');
ylabel('Magnitude');
end
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(~, ~, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global wave Fs B t
B=1;
[wave,Fs] = audioread(uigetfile({'*.mp3'})); 
t=0:(1/Fs):((size(wave)-1)/Fs);
sound(wave,Fs);
axes(handles.axes1);
N=length(wave);
a=wave(1:N);
plot(t,a);
xlabel('t/s');
ylabel('x(t)');
title('Sound Wave(TimeDomain)');
grid on;
%%用�?�显示音频信�?�的频率：
set(handles.text3,'string',sprintf('%3.0f',Fs));



function edit2_Callback(~, ~, ~)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, ~, ~)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(~, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(~, ~, ~)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, ~, ~)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(~, ~, ~)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, ~, ~)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(~, ~, ~)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, ~, ~)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(~, ~, ~)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, ~, ~)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, ~)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(~, ~, handles)%clear
    global A B C H
    A=0;
    B=0;
    C=0;
    H=0;
    set(handles.text3,'string',sprintf('%3.0f',0));%音频信�?�的抽样频率的显示也清零
try
      delete((allchild(handles.axes1)));
      axes(handles.axes1);
      grid off;
      delete((allchild(handles.axes2)));
      axes(handles.axes2);
      grid off;
      delete((allchild(handles.axes3)));
      axes(handles.axes3);
      grid off;
      delete((allchild(handles.axes4)));
      axes(handles.axes4);
      grid off;
end
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, ~, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ChosenWave wave t a f f0 A B Fs fs
A=1;%With Noise       
axes(handles.axes1);
if(ChosenWave==6)
    P0 = inputdlg({'Magnitude of Noise(dBm)'},'parameters',1,{'10'});
    P0=str2double(P0(1));
else
    P0 = inputdlg({'S/N(dB)'},'parameters',1,{'10'});
    P0=str2double(P0(1));
end

if(B==1)%mp3
 
wave=awgn(wave,10);
t=0:(1/Fs):((size(wave)-1)/Fs);
N=length(wave);
a=wave(1:N);
sound(wave,Fs);
plot(t,a);
xlabel('t/s');
ylabel('Magnitude');
grid on;
else
switch ChosenWave
    case 1
        wave=awgn(wave,P0);%%+10dB
        plot(t,wave);
        axis([0 3.5/f -1.2*a 1.2*a]);
        xlabel('t/s');
        ylabel('x(t)');
        title('With Noise TimeDomain');
        grid on;
    case 2
        wave=awgn(wave,P0);
        plot(t,wave);
        xlabel('t/s');
        ylabel('x(t)');
        title('With Noise TimeDomain');
        axis([0 3.5/f -1.2*a 1.2*a]);
        grid on;
    case 3
        wave=awgn(wave,P0);
        plot(t,wave);
        axis([0 3.5/f -1.2*a 1.2*a]);
        xlabel('t/s');
        ylabel('x(t)');
        title('With Noise TimeDomain');
        grid on;
    case 4
        wave=awgn(wave,P0);
        plot(t,wave);
        axis([0 3.5/f -1.2*a 1.2*a]);
        xlabel('t/s');
        ylabel('x(t)');
        title('Multiple Sine Signal with Noise');
        grid on;
        
    case 5
        wave=awgn(wave,P0);%S/N=30dB
        plot(t,wave);
        xlabel('t/s');
        ylabel('x(t)');
        axis auto;
        grid on;
        title('With Noise TimeDomain');
        axis ([0 3/f0 -1.2*a 1.2*a]); 
    case 6
        wave=wgn(1,1000,P0);%generates an M-by-N matrix of white Gaussian noise. Pspecifies the power of the output noise in dBW.
        t=(0:999)/fs;
        plot(t,wave);
        xlabel('t/s');
        ylabel('x(t)');
        title('With Noise TimeDomain');
        grid on;
          
end
end

guidata(hObject, handles);


% --- Executes on button press in pushbutton3.


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton3.
function pushbutton3_ButtonDownFcn(~, ~, ~)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(~, ~, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  h n
axes(handles.axes3);
stem(n,h,'.');
axis tight;
xlabel('n');ylabel('h(n)');
grid on; 
title('Window TimeDomain');

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(~, ~, handles)           %%滤波�?�的频域波形
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X Fs B fs D
D=0;
axes(handles.axes4);
     N=size(X);
     N=N(2);
 if(B==1)
     f=(0:N-1)*Fs/N;
 else
     f=(0:N-1)*fs/N;
 end
     y=fft(X); 
     mag=abs(y)*2/N;
     plot(f(1:round(N/2)),mag(1:round(N/2)));
     xlabel('Analog Freq(Hz)');
     ylabel('Magnitude');
     title('After Filter');
     grid on;


% --- Executes on button press in pushbutton3.


% --- Executes on slider movement.
function slider2_Callback(hObject, ~, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
axes(handles.axes1);
global t wave a C B
       N=length(wave);
       wave=wave(1:N);
       wmax=max(wave);
       tmax=max(t);
       c=get(hObject,'value');
       if(c==0)
           c=1;% WTF 
       end
       if(C==1)
           plot(t,wave);
           xlim([0,c*tmax]);
       grid on;
       else
       plot(t,wave);
       xlabel('t/s');
       ylabel('x(t)');
       if(B==0)
       ylim([-1.2*a,1.2*a]);
       end
       xlim([0,c*tmax]);
       grid on;
       end

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, ~, ~)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit8_Callback(hObject, ~, ~)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
global fs
str=get(hObject,'string');
fs=str2double(str);

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, ~, ~)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function slider3_Callback(hObject, ~, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
axes(handles.axes2);   %第二个图形的滑�?�
global wave fs FirstWave Fs B
       c=get(hObject,'value');
       if(c==0)
           c=1;
       end
     if(FirstWave==1)
         if(B==1)
             N=length(wave);
             f=(0:N-1)*Fs/N;
         else
             N=length(wave);
             f=(0:N-1)*fs/N;
         end
        y=fft(wave,N); 
        mag=abs(y);
        plot(f(1:round(N/2)),mag(1:round(N/2))*2/N);
        xlabel('AnalogFreq');ylabel('Magnitude');
        title('P');
        xlim([0,0.5*c*fs]);
        grid on;
     else  
        if(B==1)  
        window=boxcar(length(wave));
        N=length(wave); 
        wave=wave(1:N);
        [Pxx,f]=periodogram(wave,window,N,Fs);
        plot(f,10*log10(Pxx));
        xlim([0,0.5*c*Fs]);
        xlabel('Analog Freq (Hz)');
        ylabel('P(dB)');
        title('Soundwave P Density');
        grid on;
        else
        window=boxcar(length(wave));
        N=length(wave); 
        [Pxx,f]=periodogram(wave,window,N,fs);
        plot(f,10*log10(Pxx));
        title('P');
        grid on;
        xlim([0,0.5*c*fs]);
        xlabel('AnalogFreq(Hz)');
        ylabel('P(dB)');
        end
    end
        
% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, ~, ~)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, ~, handles)%第三个波形的滑�?�
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
       global X Fs B fs a D H wave 
       if(H==1)
           X=wave;
       end
       c=get(hObject,'value');
       if(c==0)
           c=1;
       end
       axes(handles.axes4);
       N=size(X);
       N=N(2);
if(D==1)
    if(B==1)
      t=0:(1/Fs):((N-1)/Fs);
      tmax=max(t);
      plot(t,X);
      xlim([0,c*tmax]);
      if(H==1)
          title('Denoise TimeDomain');
      else
          title('Filted TimeDomain');
      end
      grid on;
      xlabel('t/s');
      ylabel('Magnitude');
    else
      t=0:1/fs:((N-1)/fs);
      tmax=max(t);
      plot(t,X);
      title('After Filter TimeDomain');
      xlim([0,c*tmax]);
      ylim([-1.2*a,1.2*a]);
      grid on;
      xlabel('t/s');
      ylabel('Magnitude');
    end
else
    if(B==1)
     f=(0:N-1)*Fs/N;
    else
     f=(0:N-1)*fs/N;
    end
     y=fft(X); 
     mag=abs(y)*2/N;
     plot(f(1:round(N/2)),mag(1:round(N/2)));
     xlabel('Analog Freq');ylabel('Magnitude');
     title('Filted FreqDomain'); 
     xlim([0,c*0.5*fs]);
     grid on;
end
   
% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, ~, ~)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(~, ~, ~)%强制关闭键
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf);


% --- Executes on key press with focus on edit8 and none of its controls.
function edit8_KeyPressFcn(~, ~, ~)
% hObject    handle to edit8 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(~, ~, handles)%%进入图�?处�?�的按钮
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.sel_wave,'visible','off');
set(handles.text1,'visible','off');
set(handles.edit2,'visible','off');
set(handles.edit3,'visible','off');
set(handles.edit4,'visible','off');
set(handles.edit5,'visible','off');
set(handles.edit8,'visible','off');
set(handles.pushbutton1,'visible','off');
set(handles.pushbutton2,'visible','off');
set(handles.pushbutton3,'visible','off');
set(handles.pushbutton4,'visible','off');
set(handles.pushbutton5,'visible','off');
set(handles.popupmenu2,'visible','off');
set(handles.popupmenu3,'visible','off');
set(handles.pushbutton7,'visible','off');
set(handles.pushbutton8,'visible','off');
set(handles.pushbutton9,'visible','off');
set(handles.text2,'visible','off');
set(handles.text3,'visible','off');
set(handles.slider2,'visible','off');
set(handles.slider3,'visible','off');
set(handles.slider4,'visible','off');
set(handles.pushbutton15,'visible','on');
set(handles.pushbutton12,'visible','on');
set(handles.pushbutton11,'visible','off');
set(handles.pushbutton18,'visible','on');
set(handles.pushbutton19,'visible','on');
set(handles.pushbutton20,'visible','on');
set(handles.pushbutton22,'visible','on');
set(handles.pushbutton23,'visible','off');
set(handles.uipanel9,'visible','on');
set(handles.radiobutton2,'visible','on');
set(handles.radiobutton3,'visible','on');
set(handles.radiobutton4,'visible','on');
set(handles.text4,'visible','on');
set(handles.slider5,'visible','on');
set(handles.uipanel12,'visible','on');
set(handles.text6,'visible','on');
set(handles.text7,'visible','on');
set(handles.text8,'visible','on');
set(handles.text9,'visible','on');
set(handles.text10,'visible','on');
% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(~, ~, handles)   %读入一幅图�?
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I I0%  I为图形�?��?
axes(handles.axes1);%读入在第一幅图中
I=imread('haha.jpg');
zoom on;
I0=I(:,:,1);
%axes(handles.axes2)
imshow(I);
%axis('off');



% --- Executes on slider movement.
function slider5_Callback(hObject, ~, handles)   %图�?旋转函数
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
       global I;
       axes(handles.axes1);
       c=get(hObject,'value');
       if(c==0)
           c=1;
       end
       J=imrotate(I,c*360,'bilinear');
       imshow(J);
       

% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, ~, ~)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when selected object is changed in uipanel9.
% 
% function uipanel9_SelectionChangeFcn(hObject, ~, handles)
% % hObject    handle to the selected object in uipanel9 
% % eventdata  structure with the following fields (see UIBUTTONGROUP)
% %	EventName: string 'SelectionChanged' (read only)
% %	OldValue: handle of the previously selected object or empty if none was selected
% %	NewValue: handle of the currently selected object
% % handles    structure with handles and user data (see GUIDATA)
% global I
% str=get(hObject,'string');
% axes(handles.axes1);
% switch str
%     case 'Origin'
%         I=imread('haha.jpg');
%         imshow(I);
%     case 'Gray'
%         I0=I(:,:,1);%产生�?�色图�?
%         imshow(I0);
%     case '二值图�?'
%         I0=I(:,:,1);
%         I1=im2bw(I0,0.4);%二值图�?
%         imshow(I1);
%     case '索引色图�?'
%         I0=I(:,:,1);
%         x=grayslice(I0,16);
%         imshow(x,hot(16));
% end

% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(~, ~, handles)  %返回键
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.sel_wave,'visible','on');
set(handles.text1,'visible','on');
set(handles.edit2,'visible','on');
set(handles.edit3,'visible','on');
set(handles.edit4,'visible','on');
set(handles.edit5,'visible','on');
set(handles.edit8,'visible','on');
set(handles.pushbutton1,'visible','on');
set(handles.pushbutton2,'visible','on');
set(handles.pushbutton3,'visible','on');
set(handles.pushbutton4,'visible','on');
set(handles.pushbutton5,'visible','on');
set(handles.pushbutton23,'visible','on');
set(handles.popupmenu2,'visible','on');
set(handles.popupmenu3,'visible','on');
set(handles.pushbutton7,'visible','on');
set(handles.pushbutton8,'visible','on');
set(handles.pushbutton9,'visible','on');
set(handles.pushbutton18,'visible','off');
set(handles.pushbutton19,'visible','off');
set(handles.pushbutton20,'visible','off');
set(handles.pushbutton22,'visible','off');
set(handles.text2,'visible','on');
set(handles.text3,'visible','on');
set(handles.slider2,'visible','on');
set(handles.slider3,'visible','on');
set(handles.slider4,'visible','on');
set(handles.text6,'visible','off');
set(handles.pushbutton11,'visible','on');
set(handles.pushbutton12,'visible','off');
set(handles.uipanel9,'visible','off');
set(handles.radiobutton2,'visible','off');
set(handles.radiobutton3,'visible','off');
set(handles.radiobutton4,'visible','off');
set(handles.text4,'visible','off');
set(handles.slider5,'visible','off');
set(handles.uipanel12,'visible','off');
set(handles.text7,'visible','off');
set(handles.text8,'visible','off');
set(handles.text9,'visible','off');
set(handles.text10,'visible','off');
% --- Executes when selected object is changed in uipanel11.


% --- Executes when selected object is changed in uipanel11.
function uipanel11_SelectionChangeFcn(hObject, ~, handles)
% hObject    handle to the selected object in uipanel11 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global  I0
str=get(hObject,'string');
axes(handles.axes3);
switch str
    case 'poisson'
        J=imnoise(I0,'poisson');
        imshow(J);
    case 'salt & pepper'
        J=imnoise(I0,'salt & pepper', 0.02);
        imshow(J);
    case 'gaussian'
        J=imnoise(I0,'gaussian', 0.02);
        imshow(J);
end


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(~, ~, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  I
I0=I(:,:,1);
axes(handles.axes3);
J=imnoise(I0,'poisson');
imshow(J);


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(~, ~, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  I J
axes(handles.axes3);
I0=I(:,:,1);
J=imnoise(I0,'gaussian');
 imshow(J);


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(~, ~, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  I0 J
axes(handles.axes3);
J=imnoise(I0,'salt & pepper');
 imshow(J);


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(~, ~, handles)%中值滤波函数
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes4);%%J指�?�噪信�?�，K指滤波之�?�的信�?�
global J
K=medfilt2(J);
imshow(K);


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(~, ~, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I
axes(handles.axes2);
I0=I(:,:,1);
BW=edge(I0);%产生边界
imshow(BW);%�?�边缘图�?


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(~, ~, handles)%自适应滤波
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes4);
global J
K=wiener2(J,[3,3]);
imshow(K)
% --- Executes during object creation, after setting all properties.


% --- Executes during object creation, after setting all properties.
function pushbutton22_CreateFcn(~, ~, ~)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in uipanel12.
function uipanel12_SelectionChangeFcn(hObject, ~, handles)%边缘检测函数
% hObject    handle to the selected object in uipanel12 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global I
axes(handles.axes2);
str=get(hObject,'string');
switch str
    case 'sobel'
        BW=edge(rgb2gray(I),'sobel');
        imshow(BW);
    case 'prewitt'
        BW=edge(rgb2gray(I),'prewitt');
        imshow(BW);
    case 'canny'
        BW=edge(rgb2gray(I),'canny');
        imshow(BW);
end


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(~, ~, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global wave Fs B fs D f a H
H=1;
axes(handles.axes4);
D=1;
N=length(wave);
wave=wave(1:N);
lev=5;
[c,l]=wavedec(wave,lev,'sym8');
wave=wden(c,l,'minimaxi','s','sln',lev,'sym8');
if(B==1)
    t=0:(1/Fs):((size(wave)-1)/Fs);
    plot(t,wave);
    title('DeNoise TimeDomain');
    grid on;
    xlabel('t/s');
    ylabel('Magnitude');
    sound(wave,Fs);
    ylim([-1,1]);
else
    t=0:1/fs:(N-1)/fs;
    plot(t,wave);
    title('DeNoise TimeDomain');
    xlim([0,10/f]);
    ylim([-1.2*a,1.2*a]);
    grid on;
    xlabel('t/s');
    ylabel('Magnitude');
end

% --- Executes on mouse press over axes background.
function axes6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tired
axes(handles.axes6);
tired=imread('tired.jpg');
imshow(tired);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton1.
function pushbutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function axes6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on pushbutton1 and none of its controls.
function pushbutton1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
