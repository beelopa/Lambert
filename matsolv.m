function [IJK] = matsolv(PQW,TRANS)
%Transforming positions and velocity vectors found by "planetephemeris"
%function from ecliptic plane to invariable plane

      %R_IJK = R_Trans*R_PQW;
          ijk = TRANS.*PQW;
            IJK(1) = ijk(1,1)+ijk(1,2)+ijk(1,3);
            IJK(2) = ijk(2,1)+ijk(2,2)+ijk(2,3);
            IJK(3) = ijk(3,1)+ijk(3,2)+ijk(3,3);  
end

