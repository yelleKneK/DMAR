# print.dmar_tbl renders a long table with p_terms and a CI footer

    Code
      print(equivalence_smd(smd = 0.1, n_1 = 50, n_2 = 50, delta_lower = 0.5,
        delta_upper = 0.5))
    Output
       term        value 
       smd         0.1   
       t_lower     3     
       t_upper     -2    
       df          98    
       p_lower     0.0017
       p_upper     0.0242
       p_tost      0.0242
       lower_limit -0.232
       upper_limit 0.432 
       delta_lower -0.5  
       delta_upper 0.5   
       equivalent  1     
      
      Confidence level: 90%

---

    Code
      print(ci_smd(smd = 0.5, n_1 = 30, n_2 = 30))
    Output
       term        value  
       lower_limit -0.0162
       smd         0.5    
       upper_limit 1.01   
      
      Confidence level: 95%

# print.dmar_tbl renders a wide table column by column

    Code
      print(ci_mahalanobis(D2 = 103.2, n_1 = 50, n_2 = 50, p = 4))
    Output
       sample_type D2  lower_limit upper_limit F_value df_1 df_2 n_1 n_2 p
       two-sample  103 72.6        131         625     4    95   50  50  4
      
      Confidence level: 95%

# format_p renders fixed decimals with the floor label

    Code
      format_p(c(0.234567, 0.05, 0.0499, 1.010542e-05, 0, 1))
    Output
      [1] "0.2346"   "0.0500"   "0.0499"   "< 0.0001" "< 0.0001" "1.0000"  

---

    Code
      format_p(c(0.234567, 1e-04), digits_p = 2)
    Output
      [1] "0.23"   "< 0.01"

# print_anova renders an anova table with DMAR p-value formatting

    Code
      print_anova(anova(fit))
    Output
      Analysis of Variance Table
      
      Response: mpg
      
                Df    Sum Sq    Mean Sq   F value   Pr(>F)
      wt         1 847.72525 847.725250 126.04109 < 0.0001
      hp         1  83.27418  83.274183  12.38133   0.0015
      Residuals 29 195.04775   6.725785        NA     <NA>

# print_summary renders an lm summary with DMAR p-value formatting

    Code
      print_summary(summary(fit))
    Output
                    Length Class  Mode   
      call           3     -none- call   
      terms          3     terms  call   
      residuals     32     -none- numeric
      coefficients  12     -none- numeric
      aliased        3     -none- logical
      sigma          1     -none- numeric
      df             3     -none- numeric
      r.squared      1     -none- numeric
      adj.r.squared  1     -none- numeric
      fstatistic     3     -none- numeric
      cov.unscaled   9     -none- numeric

