import string

def Fib(n): #	Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(-1) = 0 and F(0) = 1.
#Fib_Seq = [ (0), 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040, 1346269 ]
    if n <= 1:
        return n #Seq in Cipher given seem to begin at the second index of fib seq and not the first which would be zero.
    else:
        return(Fib(n-1) + Fib(n-2))

def shift_by( letter, Fib_index):
    alphabet = list(string.ascii_uppercase)
    index = alphabet.index(letter)
    shifted_letter = alphabet[(index+Fib_index) % len(alphabet)] #Wraps around array
    return shifted_letter

alphabet = list(string.ascii_uppercase)


# UnknownWord = "UNICORN"
print ("Enter the word you what to encrypt:")
UnknownWord = raw_input().upper() #Get input from user
NewWordArr=[]
for i in range(len(UnknownWord)):
    NewWordArr.append(shift_by(UnknownWord[i],Fib(i+1)))

NewWord = ''.join(NewWordArr)
print ("\nThe new word is:  ")
print (NewWord)
