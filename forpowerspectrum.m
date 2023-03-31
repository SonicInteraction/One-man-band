%% magnitude spectrum
x = bass1;
wlen = 2048;            h = 0.5*wlen;
win = hanning(wlen);
pin = h;
pout = pin + h;
X = zeros(length(x),1);
while pout < length(x);
    xg = x(pin-h+1:pin+h).*win;
    XG = fft(xg);
    X(pin-h+1:pin+h) = X(pin-h+1:pin+h)+XG;
    pin = pin+h;  pout = pin+h;
end
plot(abs(XG(1:100)).^2);