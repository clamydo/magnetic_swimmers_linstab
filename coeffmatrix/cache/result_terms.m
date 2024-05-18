{(-I)*(k*Cos[\[Theta]]*Cos[\[CapitalTheta]] + k*Cos[\[Phi]]*Sin[\[Theta]]*
    Sin[\[CapitalTheta]])*SphericalHarmonicY[l, m, \[Theta], \[Phi]], 
 2*b*Cos[\[Theta]]*SphericalHarmonicY[l, m, \[Theta], \[Phi]], 
 b*m*Cos[\[Theta]]*SphericalHarmonicY[l, m, \[Theta], \[Phi]], 
 (b*Sqrt[(l - m)*(1 + l + m)]*Sin[\[Theta]]*SphericalHarmonicY[l, 1 + m, 
    \[Theta], \[Phi]])/E^(I*\[Phi]), 
 -((dt*k^2 + dr*l*(1 + l))*SphericalHarmonicY[l, m, \[Theta], \[Phi]]), 
 -(\[Kappa]*Sin[\[Theta]]^2*(-(Cos[\[Theta]]*Sin[\[Theta]]*
       (Cos[\[Theta]]*(-(Cos[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*Indexed[
                \[Sigma], {2, 3}] + Indexed[\[Sigma], {2, 1}]*Sin[
                \[CapitalTheta]]))/2 - (\[Gamma]*Cos[\[CapitalTheta]]*
            (Cos[\[CapitalTheta]]*Indexed[\[Sigma], {2, 3}] + 
             Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]]))/2) + 
        Cos[\[Phi]]*Sin[\[Theta]]*
         (-(Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*Indexed[\[Sigma], 
                {2, 3}] + Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]]))/
           2 - (\[Gamma]*Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
              Indexed[\[Sigma], {2, 3}] + Indexed[\[Sigma], {2, 1}]*
              Sin[\[CapitalTheta]]))/2))*Sin[\[Phi]]) + 
     (1 - Cos[\[Theta]]^2)*(Cos[\[Phi]]*Sin[\[Theta]]*
        ((Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}] - 
           Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*(-Indexed[\[Sigma], 
                 {1, 1}] + Indexed[\[Sigma], {3, 3}]) + 
             Indexed[\[Sigma], {3, 1}]*Sin[\[CapitalTheta]]))/2 + 
         (\[Gamma]*Cos[2*\[CapitalTheta]]*(-(Cos[\[CapitalTheta]]^2*
              Indexed[\[Sigma], {1, 3}]) + Sin[\[CapitalTheta]]*
             (Cos[\[CapitalTheta]]*(-Indexed[\[Sigma], {1, 1}] + 
                Indexed[\[Sigma], {3, 3}]) + Indexed[\[Sigma], {3, 1}]*Sin[
                \[CapitalTheta]])))/2) + (\[Gamma]*Cos[\[Theta]]*
         (Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}] + 
          Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*(Indexed[\[Sigma], {1, 
                1}] - Indexed[\[Sigma], {3, 3}]) - Indexed[\[Sigma], {3, 1}]*
             Sin[\[CapitalTheta]]))*Sin[2*\[CapitalTheta]])/2 + 
       Sin[\[Theta]]*((Cos[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
             Indexed[\[Sigma], {2, 3}] + Indexed[\[Sigma], {2, 1}]*
             Sin[\[CapitalTheta]]))/2 - (\[Gamma]*Cos[\[CapitalTheta]]*
           (Cos[\[CapitalTheta]]*Indexed[\[Sigma], {2, 3}] + 
            Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]]))/2)*
        Sin[\[Phi]]) - Cos[\[Theta]]*Cos[\[Phi]]*Sin[\[Theta]]*
      (Cos[\[Theta]]*((-(Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}]) - 
           Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*(Indexed[\[Sigma], 
                {1, 1}] - Indexed[\[Sigma], {3, 3}]) - 
             Indexed[\[Sigma], {3, 1}]*Sin[\[CapitalTheta]]))/2 + 
         (\[Gamma]*Cos[2*\[CapitalTheta]]*(-(Cos[\[CapitalTheta]]^2*
              Indexed[\[Sigma], {1, 3}]) + Sin[\[CapitalTheta]]*
             (Cos[\[CapitalTheta]]*(-Indexed[\[Sigma], {1, 1}] + 
                Indexed[\[Sigma], {3, 3}]) + Indexed[\[Sigma], {3, 1}]*Sin[
                \[CapitalTheta]])))/2) + (\[Gamma]*Cos[\[Phi]]*Sin[\[Theta]]*
         (-(Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}]) + 
          Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*(-Indexed[\[Sigma], 
                {1, 1}] + Indexed[\[Sigma], {3, 3}]) + 
            Indexed[\[Sigma], {3, 1}]*Sin[\[CapitalTheta]]))*
         Sin[2*\[CapitalTheta]])/2 + Sin[\[Theta]]*
        ((Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*Indexed[\[Sigma], 
              {2, 3}] + Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]]))/2 - 
         (\[Gamma]*Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
             Indexed[\[Sigma], {2, 3}] + Indexed[\[Sigma], {2, 1}]*
             Sin[\[CapitalTheta]]))/2)*Sin[\[Phi]]))) + 
  \[Kappa]*Cos[\[Theta]]*Cos[\[Phi]]*Sin[\[Theta]]*
   (-(Cos[\[Phi]]*Sin[\[Theta]]^2*
      (Cos[\[Theta]]*(-(Cos[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
              Indexed[\[Sigma], {2, 3}] + Indexed[\[Sigma], {2, 1}]*
              Sin[\[CapitalTheta]]))/2 - (\[Gamma]*Cos[\[CapitalTheta]]*
           (Cos[\[CapitalTheta]]*Indexed[\[Sigma], {2, 3}] + 
            Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]]))/2) + 
       Cos[\[Phi]]*Sin[\[Theta]]*
        (-(Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*Indexed[\[Sigma], {2, 
                3}] + Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]]))/2 - 
         (\[Gamma]*Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
             Indexed[\[Sigma], {2, 3}] + Indexed[\[Sigma], {2, 1}]*
             Sin[\[CapitalTheta]]))/2))*Sin[\[Phi]]) - 
    Cos[\[Theta]]*Cos[\[Phi]]*Sin[\[Theta]]*
     (Cos[\[Phi]]*Sin[\[Theta]]*
       ((Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}] - 
          Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*(-Indexed[\[Sigma], 
                {1, 1}] + Indexed[\[Sigma], {3, 3}]) + 
            Indexed[\[Sigma], {3, 1}]*Sin[\[CapitalTheta]]))/2 + 
        (\[Gamma]*Cos[2*\[CapitalTheta]]*(-(Cos[\[CapitalTheta]]^2*
             Indexed[\[Sigma], {1, 3}]) + Sin[\[CapitalTheta]]*
            (Cos[\[CapitalTheta]]*(-Indexed[\[Sigma], {1, 1}] + Indexed[
                \[Sigma], {3, 3}]) + Indexed[\[Sigma], {3, 1}]*
              Sin[\[CapitalTheta]])))/2) + 
      (\[Gamma]*Cos[\[Theta]]*(Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], 
           {1, 3}] + Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
            (Indexed[\[Sigma], {1, 1}] - Indexed[\[Sigma], {3, 3}]) - 
           Indexed[\[Sigma], {3, 1}]*Sin[\[CapitalTheta]]))*
        Sin[2*\[CapitalTheta]])/2 + Sin[\[Theta]]*
       ((Cos[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*Indexed[\[Sigma], 
             {2, 3}] + Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]]))/2 - 
        (\[Gamma]*Cos[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
            Indexed[\[Sigma], {2, 3}] + Indexed[\[Sigma], {2, 1}]*
            Sin[\[CapitalTheta]]))/2)*Sin[\[Phi]]) + 
    (1 - Cos[\[Phi]]^2*Sin[\[Theta]]^2)*
     (Cos[\[Theta]]*((-(Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}]) - 
          Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*(Indexed[\[Sigma], {1, 
                1}] - Indexed[\[Sigma], {3, 3}]) - Indexed[\[Sigma], {3, 1}]*
             Sin[\[CapitalTheta]]))/2 + (\[Gamma]*Cos[2*\[CapitalTheta]]*
          (-(Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}]) + 
           Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*(-Indexed[\[Sigma], 
                 {1, 1}] + Indexed[\[Sigma], {3, 3}]) + 
             Indexed[\[Sigma], {3, 1}]*Sin[\[CapitalTheta]])))/2) + 
      (\[Gamma]*Cos[\[Phi]]*Sin[\[Theta]]*
        (-(Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}]) + 
         Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
            (-Indexed[\[Sigma], {1, 1}] + Indexed[\[Sigma], {3, 3}]) + 
           Indexed[\[Sigma], {3, 1}]*Sin[\[CapitalTheta]]))*
        Sin[2*\[CapitalTheta]])/2 + Sin[\[Theta]]*
       ((Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*Indexed[\[Sigma], 
             {2, 3}] + Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]]))/2 - 
        (\[Gamma]*Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
            Indexed[\[Sigma], {2, 3}] + Indexed[\[Sigma], {2, 1}]*
            Sin[\[CapitalTheta]]))/2)*Sin[\[Phi]])) + 
  \[Kappa]*Cos[\[Theta]]*Sin[\[Theta]]*Sin[\[Phi]]*
   (-(Cos[\[Theta]]*Sin[\[Theta]]*Sin[\[Phi]]*
      (Cos[\[Phi]]*Sin[\[Theta]]*
        ((Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}] - 
           Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*(-Indexed[\[Sigma], 
                 {1, 1}] + Indexed[\[Sigma], {3, 3}]) + 
             Indexed[\[Sigma], {3, 1}]*Sin[\[CapitalTheta]]))/2 + 
         (\[Gamma]*Cos[2*\[CapitalTheta]]*(-(Cos[\[CapitalTheta]]^2*
              Indexed[\[Sigma], {1, 3}]) + Sin[\[CapitalTheta]]*
             (Cos[\[CapitalTheta]]*(-Indexed[\[Sigma], {1, 1}] + 
                Indexed[\[Sigma], {3, 3}]) + Indexed[\[Sigma], {3, 1}]*Sin[
                \[CapitalTheta]])))/2) + (\[Gamma]*Cos[\[Theta]]*
         (Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}] + 
          Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*(Indexed[\[Sigma], {1, 
                1}] - Indexed[\[Sigma], {3, 3}]) - Indexed[\[Sigma], {3, 1}]*
             Sin[\[CapitalTheta]]))*Sin[2*\[CapitalTheta]])/2 + 
       Sin[\[Theta]]*((Cos[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
             Indexed[\[Sigma], {2, 3}] + Indexed[\[Sigma], {2, 1}]*
             Sin[\[CapitalTheta]]))/2 - (\[Gamma]*Cos[\[CapitalTheta]]*
           (Cos[\[CapitalTheta]]*Indexed[\[Sigma], {2, 3}] + 
            Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]]))/2)*
        Sin[\[Phi]])) - Cos[\[Phi]]*Sin[\[Theta]]^2*Sin[\[Phi]]*
     (Cos[\[Theta]]*((-(Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}]) - 
          Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*(Indexed[\[Sigma], {1, 
                1}] - Indexed[\[Sigma], {3, 3}]) - Indexed[\[Sigma], {3, 1}]*
             Sin[\[CapitalTheta]]))/2 + (\[Gamma]*Cos[2*\[CapitalTheta]]*
          (-(Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}]) + 
           Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*(-Indexed[\[Sigma], 
                 {1, 1}] + Indexed[\[Sigma], {3, 3}]) + 
             Indexed[\[Sigma], {3, 1}]*Sin[\[CapitalTheta]])))/2) + 
      (\[Gamma]*Cos[\[Phi]]*Sin[\[Theta]]*
        (-(Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}]) + 
         Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
            (-Indexed[\[Sigma], {1, 1}] + Indexed[\[Sigma], {3, 3}]) + 
           Indexed[\[Sigma], {3, 1}]*Sin[\[CapitalTheta]]))*
        Sin[2*\[CapitalTheta]])/2 + Sin[\[Theta]]*
       ((Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*Indexed[\[Sigma], 
             {2, 3}] + Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]]))/2 - 
        (\[Gamma]*Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
            Indexed[\[Sigma], {2, 3}] + Indexed[\[Sigma], {2, 1}]*
            Sin[\[CapitalTheta]]))/2)*Sin[\[Phi]]) + 
    (Cos[\[Theta]]*(-(Cos[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
             Indexed[\[Sigma], {2, 3}] + Indexed[\[Sigma], {2, 1}]*
             Sin[\[CapitalTheta]]))/2 - (\[Gamma]*Cos[\[CapitalTheta]]*
          (Cos[\[CapitalTheta]]*Indexed[\[Sigma], {2, 3}] + 
           Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]]))/2) + 
      Cos[\[Phi]]*Sin[\[Theta]]*
       (-(Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*Indexed[\[Sigma], 
              {2, 3}] + Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]]))/2 - 
        (\[Gamma]*Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
            Indexed[\[Sigma], {2, 3}] + Indexed[\[Sigma], {2, 1}]*
            Sin[\[CapitalTheta]]))/2))*(1 - Sin[\[Theta]]^2*Sin[\[Phi]]^2)), 
 3*\[Gamma]*(Cos[\[Theta]]*Cos[2*\[CapitalTheta]]*Cos[\[Phi]]*Sin[\[Theta]]*
    (-(Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}]) + 
     Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*(-Indexed[\[Sigma], {1, 1}] + 
         Indexed[\[Sigma], {3, 3}]) + Indexed[\[Sigma], {3, 1}]*
        Sin[\[CapitalTheta]])) + 
   (Cos[\[Theta]]^2*(Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}] + 
      Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*(Indexed[\[Sigma], {1, 1}] - 
          Indexed[\[Sigma], {3, 3}]) - Indexed[\[Sigma], {3, 1}]*
         Sin[\[CapitalTheta]]))*Sin[2*\[CapitalTheta]])/2 + 
   (Cos[\[Phi]]^2*Sin[\[Theta]]^2*
     (-(Cos[\[CapitalTheta]]^2*Indexed[\[Sigma], {1, 3}]) + 
      Sin[\[CapitalTheta]]*(Cos[\[CapitalTheta]]*
         (-Indexed[\[Sigma], {1, 1}] + Indexed[\[Sigma], {3, 3}]) + 
        Indexed[\[Sigma], {3, 1}]*Sin[\[CapitalTheta]]))*
     Sin[2*\[CapitalTheta]])/2 - Cos[\[Theta]]*Cos[\[CapitalTheta]]*
    Sin[\[Theta]]*(Cos[\[CapitalTheta]]*Indexed[\[Sigma], {2, 3}] + 
     Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]])*Sin[\[Phi]] - 
   Cos[\[Phi]]*Sin[\[Theta]]^2*Sin[\[CapitalTheta]]*
    (Cos[\[CapitalTheta]]*Indexed[\[Sigma], {2, 3}] + 
     Indexed[\[Sigma], {2, 1}]*Sin[\[CapitalTheta]])*Sin[\[Phi]])}
