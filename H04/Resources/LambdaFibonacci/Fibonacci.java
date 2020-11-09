package Homework03;

import java.math.BigInteger;

public class Fibonacci {
    public static BigInteger fib(int n) {
        if (n == 0){
            return BigInteger.ZERO;
        } else if (n == 1){
            return BigInteger.ONE;
        }
        return fib(n - 2).add(fib(n - 1));
    }
}
