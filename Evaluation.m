function [ rmse ] = Evaluation( P )

    global Vbat0 
    Vbat = batVoltResponse( P );
    Vbat =  transpose(Vbat);

    rmse = sqrt((sum((Vbat-Vbat0).^2))./size(Vbat,1));

end

