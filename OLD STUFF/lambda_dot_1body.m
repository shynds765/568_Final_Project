function lambda_dot = lambda_dot_1body(Q1,Q2,Q3,Q4,Q5,Q6,d0,gamma,lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,pos_obj1,pos_obj2,pos_obj3,rho,u1,u2,u3,x1,x2,x3,x4,x5,x6,x_target1,x_target2,x_target3,x_target4,x_target5,x_target6)
%LAMBDA_DOT_1BODY
%    LAMBDA_DOT = LAMBDA_DOT_1BODY(Q1,Q2,Q3,Q4,Q5,Q6,D0,GAMMA,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,LAMBDA5,LAMBDA6,POS_OBJ1,POS_OBJ2,POS_OBJ3,RHO,U1,U2,U3,X1,X2,X3,X4,X5,X6,X_TARGET1,X_TARGET2,X_TARGET3,X_TARGET4,X_TARGET5,X_TARGET6)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    16-Apr-2023 15:22:48

t2 = cos(x3);
t3 = cos(x4);
t4 = cos(x5);
t5 = cos(x6);
t6 = sin(x3);
t7 = sin(x4);
t8 = sin(x5);
t9 = sin(x6);
t10 = tan(x3);
t11 = x5+x6;
t12 = x1.^2;
t13 = x1.^3;
t14 = x2.^2;
t17 = x1.*x2.*2.0;
t21 = 1.0./rho;
t22 = 1.0./x2;
t38 = 6.313477647065839e+1;
t15 = t5.^2;
t16 = t9.^2;
t18 = t5.*x2;
t19 = cos(t11);
t20 = sin(t11);
t23 = 1.0./t14;
t24 = t3.*t4;
t25 = t3.*t8;
t26 = t4.*t7;
t27 = t7.*t8;
t28 = 1.0./t6;
t29 = 1.0./t10;
t31 = -t14;
t32 = t14-1.0;
t30 = t18+1.0;
t33 = t2.*t25;
t34 = t2.*t26;
t35 = t2.*t27;
t36 = t2.*t24;
t37 = t32.*x1;
t39 = t31+1.0;
t40 = t32.^2;
t41 = -t35;
t42 = 1.0./t30;
t44 = -t36;
t45 = -t37;
t46 = sqrt(t39);
t48 = t25+t34;
t49 = t26+t33;
t43 = t42.^2;
t47 = t46.^3;
t50 = 1.0./t46;
t51 = 1.0./sqrt(t45);
t53 = t24+t41;
t54 = t27+t44;
t55 = t4.*t6.*t9.*t17.*t42;
t56 = t5.*t6.*t8.*t17.*t42;
t59 = t4.*t6.*t9.*t32.*t42;
t60 = t5.*t6.*t8.*t32.*t42;
t61 = t2.*t4.*t9.*t37.*t42;
t62 = t2.*t5.*t8.*t37.*t42;
t63 = t4.*t5.*t6.*t37.*t42;
t64 = t4.*t6.*t9.*t37.*t42;
t65 = t5.*t6.*t8.*t37.*t42;
t66 = t6.*t8.*t9.*t37.*t42;
t67 = t6.*t9.*t24.*t37.*t42;
t68 = t5.*t6.*t25.*t37.*t42;
t70 = t6.*t9.*t26.*t37.*t42;
t71 = t5.*t6.*t27.*t37.*t42;
t72 = t6.*t8.*t9.*t42.*t45;
t78 = t5.*t17.*t42.*t49;
t79 = t9.*t17.*t42.*t48;
t80 = t9.*t42.*t48.*x1.*x2.*-2.0;
t84 = t5.*t32.*t42.*t49;
t85 = t9.*t32.*t42.*t48;
t86 = t5.*t37.*t42.*t48;
t87 = t5.*t37.*t42.*t49;
t88 = t9.*t37.*t42.*t48;
t89 = t9.*t37.*t42.*t49;
t99 = t9.*t42.*t45.*t48;
t52 = t51.^3;
t57 = t5.*t22.*t47.*x1;
t69 = t4.*t5.*t6.*t9.*t37.*t43;
t73 = t6.*t8.*t15.*t37.*t43;
t74 = t4.*t6.*t16.*t37.*t43.*x2;
t75 = t6.*t8.*t9.*t18.*t37.*t43;
t76 = t4.*t5.*t6.*t9.*t43.*t45;
t77 = t6.*t8.*t15.*t43.*t45;
t81 = t5.*t17.*t42.*t53;
t82 = t9.*t17.*t42.*t54;
t83 = t9.*t42.*t54.*x1.*x2.*-2.0;
t90 = t5.*t32.*t42.*t53;
t91 = t9.*t32.*t42.*t54;
t92 = -t85;
t93 = t15.*t37.*t43.*t49;
t94 = t5.*t9.*t37.*t43.*t48;
t95 = t5.*t37.*t42.*t53;
t96 = t5.*t37.*t42.*t54;
t97 = t9.*t37.*t42.*t53;
t98 = t9.*t37.*t42.*t54;
t100 = t9.*t18.*t37.*t43.*t49;
t102 = t16.*t37.*t43.*t48.*x2;
t103 = t15.*t37.*t43.*t53;
t104 = t15.*t43.*t45.*t49;
t105 = t5.*t9.*t37.*t43.*t54;
t106 = t9.*t42.*t45.*t54;
t107 = t9.*t18.*t37.*t43.*t53;
t108 = t9.*t18.*t43.*t45.*t49;
t109 = t16.*t37.*t43.*t54.*x2;
t110 = t15.*t43.*t45.*t53;
t111 = t9.*t18.*t43.*t45.*t53;
t113 = t59+t60;
t114 = t61+t62;
t115 = t64+t65;
t121 = t63+t72;
t123 = t67+t68;
t124 = t70+t71;
t58 = -t57;
t101 = -t91;
t112 = (t9.*t12.*t32.*t38.*t52.*u1.*x2)./3.986e+4;
t116 = abs(t115);
t117 = pos_obj3+t115;
t118 = sign(t115);
t126 = t86+t97;
t127 = t89+t96;
t129 = t90+t92;
t130 = t87+t106;
t131 = t95+t99;
t146 = t55+t56+t76+t77;
t147 = t74+t75+t121;
t148 = t78+t83+t104+t105;
t149 = t80+t81+t94+t110;
t119 = abs(t117);
t120 = sign(t117);
t122 = t116.^2;
t128 = t84+t101;
t132 = abs(t130);
t133 = abs(t131);
t134 = pos_obj1+t131;
t135 = pos_obj2+t130;
t136 = sign(t130);
t137 = sign(t131);
t150 = t102+t111+t126;
t151 = t108+t109+t127;
t152 = t113.*t116.*t118.*2.0;
t153 = t114.*t116.*t118.*2.0;
t154 = t116.*t118.*t121.*2.0;
t156 = t116.*t118.*(t69+t73-t6.*t8.*t18.*t42.*x1.*2.0-t4.*t6.*t9.*t42.*x1.*x2.*2.0).*-2.0;
t162 = t116.*t118.*t147.*2.0;
t125 = t119.^2;
t138 = abs(t134);
t139 = abs(t135);
t140 = sign(t134);
t141 = sign(t135);
t142 = t132.^2;
t143 = t133.^2;
t155 = -t154;
t157 = t123.*t132.*t136.*2.0;
t158 = t124.*t133.*t137.*2.0;
t164 = -t162;
t173 = t127.*t132.*t136.*2.0;
t174 = t126.*t133.*t137.*2.0;
t175 = t128.*t132.*t136.*2.0;
t176 = t133.*t137.*(t85-t90).*-2.0;
t177 = t132.*t136.*(t88+t5.*t42.*t45.*t53).*-2.0;
t178 = t130.*t133.*t137.*2.0;
t193 = t132.*t136.*t148.*2.0;
t194 = t133.*t137.*t149.*2.0;
t195 = t132.*t136.*t151.*2.0;
t196 = t133.*t137.*t150.*2.0;
t144 = t138.^2;
t145 = t139.^2;
t159 = -t157;
t160 = t122+t142+t143;
t180 = -t178;
t200 = t152+t175+t176;
t201 = t155+t173+t174;
t202 = t156+t193+t194;
t204 = t164+t195+t196;
t161 = t125+t144+t145;
t163 = sqrt(t160);
t198 = t177+t180;
t199 = t153+t158+t159;
t165 = 1.0./t163;
t167 = t163.*x2;
t168 = sqrt(t161);
t169 = -t163;
t183 = t46.*t163.*2.0;
t188 = (lambda3.*t20.*t38.*t51.*t163.*u3)./3.986e+4;
t190 = (t19.*t29.*t38.*t51.*t163.*u3)./3.986e+4;
t191 = (lambda4.*t19.*t28.*t38.*t51.*t163.*u3)./3.986e+4;
t166 = t165.^3;
t170 = 1.0./t168;
t171 = -t168;
t179 = t37+t169;
t189 = t58+t183;
t192 = -t191;
t203 = (t165.*t200)./2.0;
t206 = (t165.*t202)./2.0;
t172 = d0+t171;
t184 = t5.*t179;
t205 = t39+t203;
t207 = -t206;
t181 = t21.*t172;
t186 = -t184;
t208 = t17+t207;
t182 = tanh(t181);
t197 = t167+t186;
t185 = t182.^2;
t187 = t185-1.0;
et7 = lambda5.*(t5.*t22.*t32.*t38.*t51.*u1.*(-1.254390366281987e-5)-(t9.*t22.*t38.*t51.*t205.*u2)./3.986e+4+(t9.*t22.*t32.*t38.*t52.*t179.*u2)./7.972e+4+(t20.*t29.*t32.*t38.*t52.*t163.*u3)./7.972e+4+(t20.*t29.*t38.*t51.*t165.*t200.*u3)./7.972e+4)-lambda2.*((t38.*t51.*u2.*(t5.*t205+t203.*x2))./3.986e+4-(t9.*t32.*t38.*t51.*u1)./7.972e+4+(t32.*t38.*t52.*t197.*u2)./7.972e+4)-lambda1.*(t112-t12.*t32.*t38.*t51.*t165.*u2.*1.505268439538384e-4-(t13.*t38.*t40.*t52.*t165.*u2)./3.986e+4+(t9.*t38.*t51.*u1.*x1.*x2)./9.965e+3+(t13.*t32.*t38.*t51.*t166.*t200.*u2)./3.986e+4);
et8 = -lambda6.*((t38.*t51.*u1.*(t5.*t22.*t47-t46.*t165.*t200))./3.986e+4+(t32.*t38.*t52.*u1.*(t57-t183))./7.972e+4-(t9.*t22.*t38.*t47.*t52.*t179.*u2)./7.972e+4-(t9.*t22.*t38.*t46.*t51.*t205.*u2)./3.986e+4)-(Q1.*(x1.*2.0-x_target1.*2.0))./2.0-(gamma.*t21.*t170.*t187.*(t113.*t119.*t120.*2.0+t128.*t139.*t141.*2.0-t138.*t140.*(t85-t90).*2.0))./4.0-(lambda3.*t19.*t32.*t38.*t52.*t163.*u3)./7.972e+4-(lambda3.*t19.*t38.*t51.*t165.*t200.*u3)./7.972e+4-(lambda4.*t20.*t28.*t32.*t38.*t52.*t163.*u3)./7.972e+4-(lambda4.*t20.*t28.*t38.*t51.*t165.*t200.*u3)./7.972e+4;
et9 = lambda2.*(t112-(t38.*t51.*u2.*(t163-t5.*t208+t206.*x2))./3.986e+4+(t9.*t38.*t51.*u1.*x1.*x2)./1.993e+4-(t38.*t52.*t197.*u2.*x1.*x2)./3.986e+4)+lambda6.*((t38.*t51.*u1.*(t50.*t167.*-2.0+t46.*t165.*t202+t5.*t46.*x1.*3.0+t5.*t23.*t47.*x1))./3.986e+4+(t9.*t38.*t50.*t51.*t179.*u2)./3.986e+4-(t38.*t52.*u1.*x1.*x2.*(t57-t183))./3.986e+4+(t9.*t23.*t38.*t46.*t51.*t179.*u2)./3.986e+4-(t9.*t22.*t38.*t46.*t51.*t208.*u2)./3.986e+4-(t9.*t38.*t46.*t52.*t179.*u2.*x1)./3.986e+4);
et10 = -lambda5.*((t5.*t23.*t38.*u1)./(t51.*3.986e+4)+(t5.*t38.*t51.*u1.*x1)./1.993e+4+(t5.*t12.*t32.*t38.*t52.*u1)./3.986e+4+(t9.*t23.*t38.*t51.*t179.*u2)./3.986e+4-(t9.*t22.*t38.*t51.*t208.*u2)./3.986e+4-(t9.*t38.*t52.*t179.*u2.*x1)./3.986e+4-(t20.*t29.*t38.*t51.*t165.*t202.*u3)./7.972e+4-(t20.*t29.*t38.*t52.*t167.*u3.*x1)./3.986e+4);
et11 = -lambda1.*((t9.*t12.*t38.*t51.*u1)./1.993e+4+(t9.*t13.*t14.*t38.*t52.*u1)./1.993e+4-(t13.*t38.*t51.*t165.*u2.*x2)./9.965e+3+(t13.*t32.*t38.*t51.*t166.*t202.*u2)./3.986e+4-(t12.^2.*t32.*t38.*t52.*t165.*u2.*x2)./1.993e+4)-(Q2.*(x2.*2.0-x_target2.*2.0))./2.0-(gamma.*t21.*t170.*t187.*(t119.*t120.*(t69+t73-t6.*t8.*t18.*t42.*x1.*2.0-t4.*t6.*t9.*t42.*x1.*x2.*2.0).*-2.0+t138.*t140.*t149.*2.0+t139.*t141.*t148.*2.0))./4.0-(lambda3.*t19.*t38.*t51.*t165.*t202.*u3)./7.972e+4-(lambda3.*t19.*t38.*t52.*t167.*u3.*x1)./3.986e+4-(lambda4.*t20.*t28.*t38.*t51.*t165.*t202.*u3)./7.972e+4;
et12 = lambda4.*t20.*t28.*t38.*t52.*t167.*u3.*x1.*(-2.508780732563974e-5);
et13 = lambda6.*((t38.*t46.*t51.*t165.*t199.*u1)./3.986e+4+(t9.*t22.*t38.*t46.*t51.*t165.*t199.*u2)./7.972e+4)-(Q3.*(x3.*2.0-x_target3.*2.0))./2.0-lambda5.*((t9.*t22.*t38.*t51.*t165.*t199.*u2)./7.972e+4-(t20.*t29.*t38.*t51.*t165.*t199.*u3)./7.972e+4+(t20.*t29.^2.*t38.*t51.*t163.*u3.*(1.0./t29.^2+1.0))./3.986e+4)-(lambda2.*t38.*t51.*u2.*((t5.*t165.*t199)./2.0+(t165.*t199.*x2)./2.0))./3.986e+4-(gamma.*t21.*t170.*t187.*(t114.*t119.*t120.*2.0+t124.*t138.*t140.*2.0-t123.*t139.*t141.*2.0))./4.0-(lambda3.*t19.*t38.*t51.*t165.*t199.*u3)./7.972e+4+(lambda4.*t2.*t20.*t28.^2.*t38.*t51.*t163.*u3)./3.986e+4;
et14 = lambda1.*t13.*t32.*t38.*t51.*t166.*t199.*u2.*(-2.508780732563974e-5)-(lambda4.*t20.*t28.*t38.*t51.*t165.*t199.*u3)./7.972e+4;
et15 = -lambda6.*((t38.*t46.*t51.*t165.*u1.*(t178+t132.*t136.*(t88+t5.*t42.*t45.*t53).*2.0))./3.986e+4+(t9.*t22.*t38.*t46.*t51.*t165.*u2.*(t178+t132.*t136.*(t88+t5.*t42.*t45.*t53).*2.0))./7.972e+4)+lambda5.*((t9.*t22.*t38.*t51.*t165.*u2.*(t178+t132.*t136.*(t88+t5.*t42.*t45.*t53).*2.0))./7.972e+4-(t20.*t29.*t38.*t51.*t165.*u3.*(t178+t132.*t136.*(t88+t5.*t42.*t45.*t53).*2.0))./7.972e+4)-(Q4.*(x4.*2.0-x_target4.*2.0))./2.0+(lambda2.*t38.*t51.*u2.*((t5.*t165.*(t178+t132.*t136.*(t88+t5.*t42.*t45.*t53).*2.0))./2.0+(t165.*x2.*(t178+t132.*t136.*(t88+t5.*t42.*t45.*t53).*2.0))./2.0))./3.986e+4;
et16 = (gamma.*t21.*t170.*t187.*(t139.*t141.*(t88+t5.*t42.*t45.*t53).*2.0+t130.*t138.*t140.*2.0))./4.0+(lambda3.*t19.*t38.*t51.*t165.*u3.*(t178+t132.*t136.*(t88+t5.*t42.*t45.*t53).*2.0))./7.972e+4+(lambda1.*t13.*t32.*t38.*t51.*t166.*u2.*(t178+t132.*t136.*(t88+t5.*t42.*t45.*t53).*2.0))./3.986e+4+(lambda4.*t20.*t28.*t38.*t51.*t165.*u3.*(t178+t132.*t136.*(t88+t5.*t42.*t45.*t53).*2.0))./7.972e+4;
et17 = t188+t192-lambda6.*((t38.*t51.*u1.*(t46.*t165.*t204-t9.*t22.*t47.*x1))./3.986e+4+(t22.*t38.*t46.*t51.*t184.*u2)./3.986e+4+(t9.*t22.*t38.*t46.*t51.*t165.*t204.*u2)./7.972e+4)-lambda1.*((t12.*t18.*t38.*t51.*u1)./1.993e+4-(t13.*t32.*t38.*t51.*t166.*t204.*u2)./3.986e+4)+lambda5.*(t190-(t9.*t22.*t38.*u1)./(t51.*3.986e+4)+(t22.*t38.*t51.*t184.*u2)./3.986e+4+(t9.*t22.*t38.*t51.*t165.*t204.*u2)./7.972e+4-(t20.*t29.*t38.*t51.*t165.*t204.*u3)./7.972e+4);
et18 = -lambda2.*((t5.*t38.*u1)./(t51.*3.986e+4)-(t38.*t51.*u2.*(-t9.*t179+(t5.*t165.*t204)./2.0+(t165.*t204.*x2)./2.0))./3.986e+4)-(Q6.*(x6.*2.0-x_target6.*2.0))./2.0+(gamma.*t21.*t170.*t187.*(t119.*t120.*t147.*-2.0+t138.*t140.*t150.*2.0+t139.*t141.*t151.*2.0))./4.0+(lambda3.*t19.*t38.*t51.*t165.*t204.*u3)./7.972e+4+(lambda4.*t20.*t28.*t38.*t51.*t165.*t204.*u3)./7.972e+4;
mt1 = [et7+et8;et9+et10+et11+et12;et13+et14;et15+et16];
mt2 = [t188+t192-lambda6.*((t38.*t46.*t51.*t165.*t201.*u1)./3.986e+4+(t9.*t22.*t38.*t46.*t51.*t165.*t201.*u2)./7.972e+4)+lambda5.*(t190+(t9.*t22.*t38.*t51.*t165.*t201.*u2)./7.972e+4-(t20.*t29.*t38.*t51.*t165.*t201.*u3)./7.972e+4)-(Q5.*(x5.*2.0-x_target5.*2.0))./2.0+(lambda2.*t38.*t51.*u2.*((t5.*t165.*t201)./2.0+(t165.*t201.*x2)./2.0))./3.986e+4+(gamma.*t21.*t170.*t187.*(t119.*t120.*t121.*-2.0+t126.*t138.*t140.*2.0+t127.*t139.*t141.*2.0))./4.0+(lambda3.*t19.*t38.*t51.*t165.*t201.*u3)./7.972e+4+(lambda1.*t13.*t32.*t38.*t51.*t166.*t201.*u2)./3.986e+4+(lambda4.*t20.*t28.*t38.*t51.*t165.*t201.*u3)./7.972e+4;et17+et18];
lambda_dot = [mt1;mt2];
