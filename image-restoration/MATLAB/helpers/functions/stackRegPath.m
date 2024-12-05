function [] = stackRegPath(s, d)
    o = loadImageStack(s);
    r = stackReg(o);
    saveImageStack(r, d);
end