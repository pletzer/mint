#include <GrExprParser.h>
#include <GrExprAdaptor.h>

void testZero() {
    
    Vec ts(11);
    ts.space(0., 1.);

    //const std::string xLineExpr = "pi/4. + 0.2*cos(5*pi/3.) - 0.2";
    const std::string xLineExpr = "0.";

    std::cout << "expr: " << xLineExpr << '\n';
    GrExprAdaptor xa(xLineExpr);
    const std::string polishPrefixExpr = xa.getPrefixExpr();
    std::cout << "Polish prefix expr: " << polishPrefixExpr << '\n';
    GrExprParser xExpr(ts.size(), polishPrefixExpr);
    xExpr.defineVariable("t", &ts);
    Vec* xs = xExpr.eval();

    for (size_t i = 0; i < ts.size(); ++i) {
        std::cout << "i = " << i << " ts = " << ts[i] << " xs = " << (*xs)[i] << '\n';
    }

}

void testPi() {
    
    Vec ts(11);
    ts.space(0., 1.);

    //const std::string xLineExpr = "pi/4. + 0.2*cos(5*pi/3.) - 0.2";
    const std::string xLineExpr = "pi";

    std::cout << "expr: " << xLineExpr << '\n';
    GrExprAdaptor xa(xLineExpr);
    const std::string polishPrefixExpr = xa.getPrefixExpr();
    std::cout << "Polish prefix expr: " << polishPrefixExpr << '\n';
    GrExprParser xExpr(ts.size(), polishPrefixExpr);
    xExpr.defineVariable("t", &ts);
    Vec* xs = xExpr.eval();

    for (size_t i = 0; i < ts.size(); ++i) {
        std::cout << "i = " << i << " ts = " << ts[i] << " xs = " << (*xs)[i] << '\n';
    }

}

void testPiDiv1() {
    
    Vec ts(11);
    ts.space(0., 1.);

    //const std::string xLineExpr = "pi/4. + 0.2*cos(5*pi/3.) - 0.2";
    const std::string xLineExpr = "pi /1";

    std::cout << "expr: " << xLineExpr << '\n';
    GrExprAdaptor xa(xLineExpr);
    const std::string polishPrefixExpr = xa.getPrefixExpr();
    std::cout << "Polish prefix expr: " << polishPrefixExpr << '\n';
    GrExprParser xExpr(ts.size(), polishPrefixExpr);
    xExpr.defineVariable("t", &ts);
    Vec* xs = xExpr.eval();

    for (size_t i = 0; i < ts.size(); ++i) {
        std::cout << "i = " << i << " ts = " << ts[i] << " xs = " << (*xs)[i] << '\n';
    }

}

void testSubtraction() {
    
    Vec ts(11);
    ts.space(0., 1.);

    const std::string xLineExpr = "pi - 3";

    std::cout << "expr: " << xLineExpr << '\n';
    GrExprAdaptor xa(xLineExpr);
    const std::string polishPrefixExpr = xa.getPrefixExpr();
    std::cout << "Polish prefix expr: " << polishPrefixExpr << '\n';
    GrExprParser xExpr(ts.size(), polishPrefixExpr);
    xExpr.defineVariable("t", &ts);
    Vec* xs = xExpr.eval();

    for (size_t i = 0; i < ts.size(); ++i) {
        std::cout << "i = " << i << " ts = " << ts[i] << " xs = " << (*xs)[i] << '\n';
    }

}


void testCos() {
    
    Vec ts(11);
    ts.space(0., 1.);

    const std::string xLineExpr = "0.2*cos(pi*t/3)";

    std::cout << "expr: " << xLineExpr << '\n';
    GrExprAdaptor xa(xLineExpr);
    const std::string polishPrefixExpr = xa.getPrefixExpr();
    std::cout << "Polish prefix expr: " << polishPrefixExpr << '\n';
    GrExprParser xExpr(ts.size(), polishPrefixExpr);
    xExpr.defineVariable("t", &ts);
    Vec* xs = xExpr.eval();

    for (size_t i = 0; i < ts.size(); ++i) {
        std::cout << "i = " << i << " ts = " << ts[i] << " xs = " << (*xs)[i] << '\n';
    }

}


void testQuarterPi() {
    
    Vec ts(11);
    ts.space(0., 1.);

    //const std::string xLineExpr = "pi/4. + 0.2*cos(5*pi/3.) - 0.2";
    const std::string xLineExpr = "pi/4.";

    std::cout << "expr: " << xLineExpr << '\n';
    GrExprAdaptor xa(xLineExpr);
    const std::string polishPrefixExpr = xa.getPrefixExpr();
    std::cout << "Polish prefix expr: " << polishPrefixExpr << '\n';
    GrExprParser xExpr(ts.size(), polishPrefixExpr);
    xExpr.defineVariable("t", &ts);
    Vec* xs = xExpr.eval();

    for (size_t i = 0; i < ts.size(); ++i) {
        std::cout << "i = " << i << " ts = " << ts[i] << " xs = " << (*xs)[i] << '\n';
    }

}

void testFull() {
    
    Vec ts(11);
    ts.space(0., 1.);

    //const std::string xLineExpr = "pi/4. + 0.2*cos(5*pi*t/3.) - 0.2";
    const std::string xLineExpr = "0.2*cos(5*pi*t/3.) - 0.2";

    std::cout << "expr: " << xLineExpr << '\n';
    GrExprAdaptor xa(xLineExpr);
    const std::string polishPrefixExpr = xa.getPrefixExpr();
    std::cout << "Polish prefix expr: " << polishPrefixExpr << '\n';
    GrExprParser xExpr(ts.size(), polishPrefixExpr);
    xExpr.defineVariable("t", &ts);
    Vec* xs = xExpr.eval();

    for (size_t i = 0; i < ts.size(); ++i) {
        std::cout << "i = " << i << " ts = " << ts[i] << " xs = " << (*xs)[i] << '\n';
    }

}



int main(int argc, char** argv) {

    testZero();
    testPi();
    testPiDiv1();
    testSubtraction();
    testCos();
    testQuarterPi();
    testFull();

    return 0;
}
