/************************************************************
  一次方程/计算器
  copyright (C) 2009-2020 by Huang Yuanbing
  version 2.00
  mail to: bailuzhou@163.com
  free use for non-commercial purposes

  some bug for input large number
 *************************************************************/

# include <limits.h>
# include <deque>
# include <limits>
# include <cctype>
# include <cstdio>
# include <memory>
# include <cstring>
# include <cmath>
# include <cstdlib>
using namespace std;

# define X  'x'
# define MAXC 256
# define ISOPER(c) (c == '+' || c == '-' || c == '*' || c == '/' || c == '^')

#if _WIN32
typedef __int64 int64;
#else
typedef long long int64;
#endif

typedef int64 stype;

bool showLog = false;
bool showDot = false;
bool isLinear = true;

stype gcd(stype a, stype b)
{
	if (b == 0)
		return a;
	return gcd(b, a % b);
}


/*
   analzye the expression is valid or not
   and replaceChar some char with mathematical symbol
   */

class Analyzer
{
private :
	bool validCharSet[MAXC];
	int charNums[MAXC];
	char inputExpr[MAXC * 2];
	char convertExpr[MAXC * 2];
	bool isValid;

	void initCharSet( )
	{
		const char* validChar = "xy0123456789+-*/^()[]{} 	.%=\r\n";
		memset(validCharSet, false, sizeof(validCharSet));
		memset(charNums, 0, sizeof(charNums));
		memset(convertExpr, 0, sizeof(convertExpr));
		for (int i = 0; validChar[i]; i++)
			validCharSet[validChar[i]] = true;
	}

public :

	char* getExp( )
	{
		return convertExpr;
	}

	bool isEquation( )
	{
		return charNums['='] == 1 && convertExpr[0] != '=';
	}

	Analyzer(const char* epstr)
	{
		isValid = true;
		initCharSet( );
		if (epstr) {
			int leng = strlen(epstr);
			if (leng > MAXC) {
				printf("Warn: too many chars\n");
				leng = MAXC - 1;
			}

			for (int i = 0; i < leng; i++) {
				char c = epstr[i];
				if (!validCharSet[c]) {
					isValid = false;
					printf("Error: invalid char %c\n", c);
					break;
				} else {
					inputExpr[i] = c;
					charNums[c] ++;
				}
			}
			if (isValid)
				inputExpr[leng] = '\0';
		}
	}

	void removeSpace( )
	{
		/* 		if (!strchr(convertExpr, ' ') &&
				!strchr(convertExpr, '\t')) {
				return;
				} */

		char lastc = convertExpr[0];
		int nchars = isspace(lastc) ? 0 : 1;

		for (int i = 1; convertExpr[i]; i++) {
			char currc = convertExpr[i];
			if (!isspace(currc))
				convertExpr[nchars++] = currc;
			lastc = currc;
		}

		if (nchars == 0) {
			isValid = false;
		} else {
			convertExpr[nchars] = '\0';
		}

		if (showLog) {
			printf("in removeSpace line %d \n", __LINE__);
			puts(convertExpr);
		}
	}

	void replaceSpace( )
	{
		if (!strchr(inputExpr, ' ') && !strchr(inputExpr, '\t')) {
			return;
		}
		//&& !strchr(inputExpr, '\n') && !strchr(inputExpr, '\r')

		char lastc = inputExpr[0];
		int nchars = isspace(lastc) ? 0 : 1;

		//remove continous space and keep only one
		for (int i = 1; inputExpr[i]; i++) {
			char currc = inputExpr[i];
			if (!isspace(currc) || !isspace(lastc))
				inputExpr[nchars++] = currc;
			lastc = currc;
		}

		//remove tail space
		while (nchars > 0 && isspace(inputExpr[nchars - 1]))
			nchars--;

		if (nchars == 0 || inputExpr[0] == '\0') {
			isValid = false;
		} else {
			inputExpr[nchars] = '\0';
			//bug fix 2 2
			lastc = inputExpr[0];
			for (int i = 0; i < inputExpr[i]; i++) {
				char currc = inputExpr[i];
				if (isspace(currc) && isdigit(lastc) && isdigit(inputExpr[i + 1])) {
					isValid = false;
					puts("Error: two number is splited by space");
					break;
				}
				lastc = currc;
			}
		}

		if (showLog) {
			printf("in replaceSpace line %d \n", __LINE__);
			puts(inputExpr);
		}
	}

	void replaceBracket( )
	{
		if (showLog) {
			printf("in replaceBracket line %d \n", __LINE__);
			puts(inputExpr);
		}
		if (charNums['('] != charNums[')']) {
			puts("Error: bracket '(' not equal to ')' ");
			isValid = false;
		} else if (charNums['='] > 1) {
			puts("Error: too many =");
			isValid = false;
		} else if (inputExpr[0] == '=') {
			puts("Error: start with = ");
			isValid = false;
		}
		//remove empty and table
		int brackets = 0, nchars = 0;
		for (int i = 0; i < inputExpr[i]; i++) {
			char c = inputExpr[i];
			if (!isspace(c))
				convertExpr[nchars++] = c;
			if (c == '(')
				brackets ++;
			else if (c == ')') {
				brackets--;
				if (brackets < 0) {
					isValid = false;
					puts("Error: too many ) before (");
					break;
				}
			}
		}

		if (nchars == 0) {
			puts("Warn: only empty space");
			isValid = false;
		}

		replaceChar('[', '(');
		replaceChar(']', ')');
		replaceChar('{', '(');
		replaceChar('}', ')');
		convertExpr[nchars] = 0;
	}

	void replaceChar(char s, char d)
	{
		char *ps = strchr(convertExpr, s);
		while (ps) {
			*ps	= d;
			ps = strchr(convertExpr, s);
		}
	}

	//0.2 + 1.112%  1/2.5 1/2.5/2  16.9(2x+1)
	//replaceChar x.y with (x + y/10), ex 1.2 --> (12/10)
	void replaceDot( )
	{
		if (showLog) {
			printf("in replaceDot line %d \n", __LINE__);
			puts(convertExpr);
		}

		char* pd = strchr(convertExpr,'.');
		while (pd) {
			if (isdigit(pd[-1]) && isdigit(pd[1])) {

				stype n1 = 0, n2 = 0, m1 = 1;
				char* pe = pd + 1;
				for (; isdigit(*pe); pe++) {
					n1 = n1 * 10 + *pe - '0';
					m1 *= 10;
				}

				char* ps = pd - 1;
				while (isdigit (*ps))
					ps--;

				for (char* pc = ps + 1; pc < pd; pc++)
					n2 = n2 * 10 + *pc - '0';

				n1 = m1 * n2 + n1;
				char pn[100] = {0};
				if (*pe != '(')
					sprintf(pn, "(%lld/%lld)", n1, m1);
				else
					sprintf(pn, "(%lld/%lld)*", n1, m1);
				// ps pd    pe
				//x + 2.1 - x = 0
				int leng = strlen(pn) + 0;
				memmove(ps + leng + 1, pe, strlen(pe) + 3);
				memcpy(ps + 1, pn, leng);
			} else {
				isValid = false;
				puts("Error: not number near .");
				break;
			}
			pd = strchr(pd + 1 ,'.');
		}
	}

	//3.1% , x%, 20%  1/25%, (20+1)%
	//replaceChar x% with (x/100), ex 20% --> (20 / 100)
	void replacePercent( )
	{
		if (showLog) {
			printf("in replacePercent line %d \n", __LINE__);
			puts(convertExpr);
		}

		char* pp = strchr(convertExpr,'%');

		while (pp) {
			if ( isdigit(pp[-1]) ) {
				char* ps = pp - 1;
				char* pe = pp + 1;
				while (isdigit(*ps))
					ps--;
				char pn[MAXC] = {0};
				stype dn = atoi(ps + 1);
				sprintf(pn, "(%lld/100)", dn);
				int leng = strlen(pn);
				memmove(ps + leng + 1, pe, strlen(pe) + 3);
				memcpy(ps + 1, pn, leng);
			} else if (pp[-1] == X || pp[-1] == ')') {
				char* ps = pp - 1;
				char* pe = pp + 1;
				char pn[MAXC] = {0};
				if (pp[-1] == ')') {
					while (*ps != '(')
						ps--;
				}
				pp[0] = 0;
				sprintf(pn, "(%s/100)", ps);
				int leng = strlen(pn);
				memmove(ps + leng, pe, strlen(pe) + 3);
				memcpy(ps, pn, leng);
			} else {
				isValid = false;
				puts("Error: invalid operator before %");
				break;
			}

			pp = strchr(pp + 3,'%');
		}
	}

	//the unitary  (operator x with be replaced with (0 operator x)
	// (-2) --> (0 - 2)
	bool addStarZero( )
	{
		if (showLog) {
			printf("in addStarZero line %d \n", __LINE__);
			puts(convertExpr);
		}

		char buf[MAXC + MAXC] = {0};
		char* pc = buf;

		for (int i = 0; convertExpr[i]; i++) {
			char c = convertExpr[i];
			//bug fix : x = 2(x+1) or x = 2x-1
			if (isdigit(pc[-1]) && (c == X || c == '('))
				*pc++ = '*';
			//bug fix : 2(x) (2+3)x x(2+1)
			else if (c == X && pc[-1] == ')')
				*pc++ = '*';
			//bug fix : x = -1 or x = -(-x) or -x = 2
			else if ((c == '+' || c == '-') &&
					(pc[-1] == '(' || pc[-1] == '=' || 0 == buf[0])) {
				*pc++ = '0';
			}
			*pc++ = c;
		}

		strcpy(convertExpr, buf);

		return true;
	}

	int getComplexty(char c)
	{
		int complexty = -12345678;

		if (c == '(')
			complexty = 2;
		else if (c == ')')
			complexty = -2;
		else if (isdigit(c) || c == X)
			complexty = 1;
		else if (ISOPER(c))
			complexty = -1;
		else if (c == '=')
			complexty = 0;
		return complexty;
	}

	//digit, X, oper, ( )
	bool checkComplexty( )
	{
		if (showLog) {
			printf("in checkComplexty line %d \n", __LINE__);
			puts(convertExpr);
		}

		char lastc = convertExpr[0];
		int complx = getComplexty(lastc);

		for (int i = 1; convertExpr[i] && isValid; i++) {
			char currc = convertExpr[i];
			bool currd = isdigit(currc), lastd = isdigit(lastc);
			bool curro = ISOPER(currc), lasto = ISOPER(lastc);

			if ((curro || currc == '=') && convertExpr[i + 1] == 0) {
				printf("Error: wrong %c in last\n", currc);
				isValid = false;
			} else if (lastc == '(') {
				if (curro || currc == ')') {
					printf("Error: wrong %c before %c\n", lastc, currc);
					isValid = false;
				}/* else if (isspace(currc) && ')' == convertExpr[i + 1]) {
					isValid = false;
					puts("Error: ()");
					}*/
			} else if (lastc == ')') {
				if (currd || currc == '(') {
					printf("Error: wrong %c before %c\n", lastc, currc);
					isValid = false;
				}
			} else if (lastc == X) {
				if (currc == X || currd || currc == '(') {
					printf("Error: wrong %c before %c\n", lastc, currc);
					isValid = false;
				}
			} else if (lastd) {
				if (currc == X ) {
					printf("Error: wrong %c before %c\n", lastc, currc);
					isValid = false;
				}/* else if (isspace(currc) && isdigit(convertExpr[i + 1])) {
					isValid = false;
					puts("Error: two digit split by space");
					}*/
			} else if (lasto || lastc == '=') {
				if (lastc == '=' && (complx < 0 || complx > 1)) {
					puts("Error: complexty failed before = ");
					isValid = false;
				}
				if (curro || currc == ')' || currc == '=') {
					printf("Error: wrong %c before %c\n", lastc, currc);
					isValid = false;
				}
			} else {
				printf("invlaid input %c or not precommand\n", currc);
			}

			if (!(currd && lastd)) {
				complx += getComplexty(currc);
				if (complx < 0) {
					isValid = false;
					printf("Error: not invlaid express for complexty = %d\n", complx);
				}
			}
			lastc = currc;
		}

		return isValid;
	}

	bool getState( )
	{
		replaceSpace( );
		if (isValid)
			replaceBracket();
		if (isValid)
			replaceDot();
		if (isValid)
			replacePercent();
		if (isValid)
			addStarZero();
		if (isValid)
			checkComplexty();

		return isValid;
	}
};

/**
  class Fraction: numerator / denominator
  and Fraction can support arithmetic + - / *
  */
class Fraction
{
private:
	stype numerator;
	stype denominator;

public:

	Fraction( )
	{
		numerator = 0;
		denominator = 1;
	}

	Fraction(stype n, stype m = 1)
	{
		numerator = n;
		denominator = m;
		if (m == 0) {
			puts("the denominator can not be 0");
			exit(2);
		} else if (m < 0) {
			n = -n;
			m = -m;
		}
		if (m > 1)
			Reduction();
	}

	Fraction& operator = (stype m)
	{
		numerator = m;
		denominator = 1;
		return *this;
	}

	void setReverse( )
	{
		numerator = -numerator;
	}

	bool operator == (stype mn) const
	{
		return numerator == mn * denominator;
	}

	bool operator != (stype mn) const
	{
		return numerator != mn * denominator;
	}

	bool operator > (stype mn) const
	{
		return numerator > mn * denominator;
	}

	// numerator/denominator + r.numerator/r.denominator
	Fraction& operator += (const Fraction &r)
	{
		//overflow
		stype factor = gcd(denominator, r.denominator);
		if (factor > 1) {
			numerator = r.numerator * (denominator / factor) +
				(r.denominator / factor) * numerator;
			denominator *= r.denominator / factor;
		} else {
			numerator = r.numerator * denominator + r.denominator * numerator;
			denominator *= r.denominator;
		}
		Reduction();
		return *this;
	}

	Fraction operator + (const Fraction &r)
	{
		Fraction tmp = *this;
		return tmp += r;
	}

	// numerator/denominator - r.numerator/r.denominator
	Fraction& operator -= (const Fraction &r)
	{
		stype factor = gcd(denominator, r.denominator);
		if (factor > 1) {
			numerator = -r.numerator * (denominator / factor) +
				(r.denominator / factor) * numerator;
			denominator *= r.denominator / factor;
		} else {
			numerator = -r.numerator * denominator + r.denominator * numerator;
			denominator *= r.denominator;
		}
		Reduction();
		return *this;
	}

	Fraction operator - (const Fraction &r)
	{
		Fraction tmp = *this;
		return tmp -= r;
	}

	Fraction& operator *= (const Fraction &r)
	{
		stype factor1 = gcd(denominator, r.numerator);

		stype factor2 = gcd(numerator, r.denominator);

		if (factor1 == 0 || factor2 == 0) {
			*this = 0;
		} else {
			denominator /= factor1;
			numerator /= factor2;
			denominator *= r.denominator / factor2;
			numerator *= r.numerator / factor1;
			Reduction();
		}
		return *this;
	}

	Fraction operator * (const Fraction &r)
	{
		Fraction tmp = *this;
		return tmp *= r;
	}

	Fraction operator ^ (const Fraction &r)
	{
		Fraction tmp = 1;
		for (int i = 0; i < r.numerator; i++)
			tmp *= *this;
		return tmp;
	}

	Fraction& operator /= (const Fraction &r)
	{
		if (r != 0) {
			stype factor1 = gcd(denominator, r.denominator);
			stype factor2 = gcd(numerator, r.numerator);

			if (factor1 == 0 || factor2 == 0) {
				*this = 0;
			} else {
				denominator /= factor1;
				numerator *= r.denominator / factor1;

				numerator /= factor2;
				denominator *= r.numerator / factor2;
				Reduction();
			}
		} else {
			printf("%lld / %lld\n", r.denominator, r.numerator);
			puts("devide 0!!!");
			exit(1);
		}

		return *this;
	}

	Fraction operator / (const Fraction &r) const
	{
		Fraction tmp = *this;
		return tmp /= r;
	}

	void Reduction( )
	{
		stype factor = gcd(denominator, numerator);
		if (factor != 1) {
			denominator /= factor;
			numerator /= factor;
		}
		if (denominator < 0) {
			denominator = -denominator;
			numerator = -numerator;
		}
	}

	void print(bool newline = true) const
	{
		if (denominator == 1)
			printf("%lld", (int64)numerator);
		else {
			if (!showDot)
				printf("%lld/%lld", (int64)numerator, (int64)denominator);
			else
				printf("%.3lf", 1.0*numerator / denominator);
		}
		if (newline)
			putchar('\n');
	}
};

/**
  class Axb save the linear express ax + b
  a and b all are Fraction
  and Axb can support arithmetic + - / *
  */
class Axb
{
private :
	Fraction a, b;

public :

	Axb(stype m, stype n)
	{
		a = m;
		b = n;
	}

	Axb& operator += (const Axb &r)
	{
		a += r.a;
		b += r.b;
		return *this;
	}

	Axb operator + (const Axb &r)
	{
		Axb tmp = *this;
		return tmp += r;
	}

	Axb& operator -= (const Axb &r)
	{
		a -= r.a;
		b -= r.b;
		return *this;
	}

	Axb operator - (const Axb &r)
	{
		Axb tmp = *this;
		return tmp -= r;
	}
	//r must be number with a = 0
	Axb& operator *= (const Axb &r)
	{
		if (a != 0 && r.a != 0) {
			printf("Warn : not linear equation for ");
			print(false);
			printf(" * ");
			r.print(true);
			isLinear = false;
		} else if (r.a == 0) {
			a *= r.b;
			b *= r.b;
		} else {
			a = b * r.a;
			b *= r.b;
		}

		return *this;
	}

	Axb operator * (const Axb &r)
	{
		Axb tmp = *this;
		return tmp *= r;
	}
	//r must be number with a = 0 and b != 0
	Axb& operator /= (const Axb &r)
	{
		if (r.a != 0) {
			printf("Warn : not linear equation for ");
			print(false);
			printf(" / ");
			r.print(true);
			isLinear = false;
		} else {
			a /= r.b;
			b /= r.b;
		}
		return *this;
	}

	Axb operator / (const Axb &r)
	{
		Axb tmp = *this;
		return tmp /= r;
	}

	Axb operator ^ (const Axb &r)
	{
		Axb tmp = *this;
		return tmp ^= r;
	}

	Axb operator ^= (const Axb &r)
	{
		if (r.a != 0 || a != 0) {
			printf("Warn : not linear equation for ");
			print(false);
			printf(" ^ ");
			r.print(true);
			isLinear = false;
		} else {
			a = a ^ r.b;
			b = b ^ r.b;
		}
		return *this;
	}

	void solveEquation( )
	{
		if (a == 0) {
			if (b != 0)
				puts("方程没有解");
			else
				puts("方程无穷多解");
		} else {
			Fraction x = b / a;
			printf("x = ");
			x.setReverse( );
			x.print(true);
		}
	}

	void print(bool newline = false) const
	{
		if (a != 0) {
			if (a != 1 && a != -1)
				a.print(false);
			else if (a == -1)
				putchar('-');
			putchar(X);
			if (b > 0) {
				putchar('+');
				b.print(false);
			} else if (b != 0)
				b.print(false);
		} else {
			b.print(false);
		}
		if (newline) {
			putchar('\n');
		}
	}
};

struct Oper
{
	enum OP
	{
		OADDSUB = 1,
		OMULDIV = 2,
		OLEFTBR = 3,
		OREFTBR = 0
	};

	char Operator;
	int priority;

	Oper(char O)
	{
		Operator = O;
		if (Operator == '+' || Operator == '-')
			priority = OADDSUB;
		else if (Operator == '*' || Operator == '/' || Operator == '^')
			priority = OMULDIV;
		else if (Operator == '(')
			priority = OLEFTBR;
		else if (Operator == ')')
			priority = OREFTBR;
	}
};

/**
  calculate the input express
  use two stack to simulate
  */
class Calculator
{
private :
	char expre[MAXC + MAXC];
	deque <Axb> daxb;
	deque <Oper> dope;

public :

	Calculator(char* pexpre)
	{
		if (showLog) {
			printf("Calculator : ");
			puts(pexpre);
		}
		memcpy(expre, pexpre, strlen(pexpre) + 1);
		removeEqual();
	}

	bool removeEqual( )
	{
		char* pe = strchr(expre, '=');
		/**if (pe == expre) {
		  puts("Warn  : with only = ");
		  return false;
		 }*/
		if (pe) {
			int leng = strlen(pe);
			memmove(pe + 2, pe + 1, leng + 2);
			pe[0] = '-';
			pe[1] = '(';
			pe[leng + 1] = ')';
			pe[leng + 2] = 0;
		}
		if (showLog) {
			printf("Calculator removeEqual : ");
			puts(expre);
		}
		return true;
	}

	bool adjustQueue( )
	{
		if (daxb.size() >= 2 && dope.size() < 1) {
			puts("stack wrong");
			return false;
		}

		if (daxb.size() < 2 || dope.size() < 1)
			return false;

		Axb& a = daxb.back();
		daxb.pop_back();
		Axb& b = daxb.back();
		daxb.pop_back();
		char oper = dope.back().Operator;

		if (showLog) {
			printf("a = ");
			a.print();
			printf(" %c b = ", oper);
			b.print(true);
		}

		if (oper == '*')
			daxb.push_back(a * b);
		else if (oper == '/')
			daxb.push_back(b / a);
		else if (oper == '+')
			daxb.push_back(a + b);
		else if (oper == '-')
			daxb.push_back(b - a);
		else if (oper == '^')
			daxb.push_back(b ^ a);

		if (showLog) {
			printf("new = ");
			daxb.back().print(true);
		}
		dope.pop_back();

		return true;
	}

	stype atoint64(char* &exp)
	{
		stype num = 0;
		while (isdigit(*exp)) {
			if (num > LONG_LONG_MAX / 10 || num < LONG_LONG_MIN / 10) {
				puts("Warn: atoint64 number overflow");
				break;
			}
			num = 10 * num + *exp++ - '0';
		}
		exp--;
		return num;
	}

	Axb expval( )
	{
		char* exp = expre;
		if (showLog){
			puts(exp);
		}

		while (*exp) {
			char c = *exp;
			if (isdigit(c) || c == X) {
				if (c != X)
					daxb.push_back(Axb(0, atoint64(exp)));
				else
					daxb.push_back(Axb(1, 0));
				if (daxb.size() >= 2 && dope.back().priority == dope.back().OMULDIV)
					adjustQueue();
			} else if (ISOPER(c)) {
				dope.push_back(Oper(c));
				if (daxb.size() >= 2 &&
						dope[dope.size() - 2].priority != dope.back().OLEFTBR &&
						dope.back().priority <= dope[dope.size() - 2].priority) {
					dope.pop_back();
					adjustQueue();
					dope.push_back(Oper(c));
				}
			} else if (c == '(' || c == ')') {
				if (c == ')') {
					while (dope.back().Operator != '(')
						adjustQueue();
					dope.pop_back();
				} else
					dope.push_back(Oper(c));
				while (dope.size() > 0 && dope.back().priority == dope.back().OMULDIV)
					adjustQueue();
			}
			exp++;
		}

		//	if (dint.size() == 1)
		//		adjustQueue();

		while (daxb.size() >= 2)
			adjustQueue();
		return daxb[0];
	}
};

static void printInfo( )
{
	puts("---------------------------------------------------------------");
	puts("1.求解带x的一元一次线性方程: 2x/3 = (x - 2) * 3.5 - 10");
	puts("2.四则运算计算器: (1.2 + 2.8)/3.1 - 25% * 6");
	puts("3.输入 dot 显示小数结果");
	puts("4.输入 log 显示计算过程");
	puts("5.Copyright by Huang Yuanbing 2009 - 2018 22738078@qq.com");

#ifdef _MSC_VER
	printf("Compiled by MS/vc++ %d", _MSC_VER);
#elif __GNUC__
	printf("Compiled by GNU/g++ %d.%d.%d",
			__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif

#if _M_AMD64 || __x86_64__
	printf(" on Windows 64 bit");
#endif
	printf(" on %s %s\n", __TIME__, __DATE__);
	puts("---------------------------------------------------------------\n");
}

void unitTest()
{
	Axb a(1, 1);
	int m, n;
	while (scanf ("%d %d", &m, &n) == 2) {
		Axb t(m, n);
		t.print();

		a += t;
		a.print();

		a -= t;
		a.print();

		a /= Axb(0,m);
		a.print();

		a *= Axb(0,m);
		a.print();
	}
}

int main( )
{
	printInfo( );

	int ncase = 1;
	char caxb[MAXC] = {0};
	while ( fgets(caxb, MAXC, stdin) ) {
		caxb[strlen(caxb) - 1] = 0;
		if (strcmp(caxb, "log") == 0) {
			if (showLog)
				puts("disable log");
			else
				puts("enable log");
			showLog = !showLog;
			continue;
		} else if (strcmp(caxb, "dot") == 0) {
			showDot = !showDot;
			if (showDot)
				puts("display format a.b");
			else
				puts("display format a/b");
			continue;
		} else if (strcmp(caxb, "exit") == 0) {
			return 0;
		}

		isLinear = true;
		Analyzer a(caxb);
		if ( a.getState() ) {
			printf("case %d result : ", ncase++);
			Calculator c(a.getExp());
			Axb ab = c.expval();
			if (!isLinear) {
				puts("not linear equation");
			} else if ( !a.isEquation() )
				ab.print();
			else
				ab.solveEquation();
			putchar('\n');
		}

		putchar('\n');
		memset(caxb, 0, sizeof(caxb));
	}

	return 0;
}

// x / 0.5
