function err=errofdata(material)
    Targert=readmatrix('Moment_abaqus.csv');
    PRED=myNeuralNetworkFunction(material);
    err1=0.1*sum((Targert(1:8)-PRED(1:8)).^2);
    err2=1*sum((Targert(9:20)-PRED(9:20)).^2);
    err3=0.01*sum((Targert(21:31)-PRED(21:31)).^2);
    err=sqrt((err1+err2+err3)/31);
end