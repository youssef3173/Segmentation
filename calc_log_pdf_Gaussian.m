function Lp = calc_log_pdf_Gaussian( I, mu, sigma )
    d = 10;
    s = size( I);
    Lp = zeros( s(1), s(2));
    for x = 1:s(1)
        for y = 1:s(1)
            
            vect = zeros( 1, 10);
            for i = 1:10
               vect( i) = I(x,y,i); 
            end
            Lp( x, y) = -log( (2*pi)^(d/2)*det(sigma) ) -0.5*(vect - mu)*pinv( sigma)*(vect - mu)';
        end
    end
end