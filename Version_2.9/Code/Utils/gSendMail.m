function gSendMail(varargin)
% Define these variables appropriately:
mail = 'gilimexico@gmail.com'; %Your GMail email address
password = 'mexicocuba'; %Your GMail password

% Then this code will set up the preferences properly:
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

% sendmail('giladliberman@gmail.com','Test subject','Test message',{});
sendmail(varargin{:});

% %%
% sendmail('giladliberman@gmail.com','Test subject','Test message',{'Simple_4.jpg'});
% %%
% D=dir('*.jpg');
% FNs={D.name};
% %%
% sendmail('giladliberman@gmail.com','Test subject','Test
% message',FNs(~strhas(FNs,'simple'))');