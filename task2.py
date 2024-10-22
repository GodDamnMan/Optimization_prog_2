
import numpy
from numpy.linalg import norm
from task1 import Simplex


class IteriorPoint:

    def __init__(self, coefs, constraints, right_hand_side, epsilon):
        self.solution_l1 = None
        self.solution_l2 = None
        self.epsilon = epsilon
        self.right_hand_side = right_hand_side
        self.constraints = constraints
        self.coefs = coefs
        self.alpha1 = 0.5
        self.alpha2 = 0.9
        self.iter_val = 1

    def find_initial_solution(self):
        initial_solution = []
        for item in self.coefs:
            if item != 0:
                initial_solution.append(1)

        for i in range(len(self.constraints)):
            row = self.constraints[i]
            res = self.right_hand_side[i]
            for item in row:
                if(item != 1):
                    res -= item
            initial_solution.append(res)
        
        self.solution_l1 = initial_solution
        self.solution_l2 = initial_solution    

    def iterior_point_method(self):
        self.find_initial_solution()
        #TODO понять когда у нас метод не будет работать
        #TODO разобраться что и как нужно выводить у симплекса
        while True:

            #calculating previous solution and diagonal matrix for lambda1 and lambda2
            prev_sol_l1 = self.solution_l1
            diagonal_l1 = numpy.diag(self.solution_l1)    

            prev_sol_l2 = self.solution_l2
            diagonal_l2 = numpy.diag(self.solution_l2)

            #calculating matrices for lambda1 and lambda2
            AA_l1 = numpy.dot(self.constraints,diagonal_l1)
            cc_l1 = numpy.dot(diagonal_l1, self.coefs)

            AA_l2 = numpy.dot(self.constraints,diagonal_l2)
            cc_l2 = numpy.dot(diagonal_l2, self.coefs)


            I = numpy.eye(len(self.coefs))


            #components to calculate P for lambda1 and lambda2
            F_l1 = numpy.dot(AA_l1, numpy.transpose(AA_l1))
            FI_l1 = numpy.linalg.inv(F_l1)
            H_l1 = numpy.dot(numpy.transpose(AA_l1), FI_l1)
            
            
            F_l2 = numpy.dot(AA_l2, numpy.transpose(AA_l2))
            FI_l2 = numpy.linalg.inv(F_l2)
            H_l2 = numpy.dot(numpy.transpose(AA_l2), FI_l2)

            #calculating P for lambda1 and lambda2
            P_l1 = numpy.subtract(I, numpy.dot(H_l1, AA_l1))
            P_l2 = numpy.subtract(I, numpy.dot(H_l2, AA_l2))

            #calculating cp for lambda1 and lambda2
            cp_l1 = numpy.dot (P_l1, cc_l1)
            cp_l2 = numpy.dot (P_l2, cc_l2)


            
            nu_l1 = numpy.absolute(numpy.min(cp_l1))
            y_l1 = numpy.add(numpy.ones(len(self.coefs), float), (self.alpha1/nu_l1) * cp_l1)

            nu_l2 = numpy.absolute(numpy.min(cp_l2))
            y_l2 = numpy.add(numpy.ones(len(self.coefs), float), (self.alpha2/nu_l2) * cp_l2)
            


            yy_l1 = numpy.dot(diagonal_l1, y_l1)
            yy_l2 = numpy.dot(diagonal_l2, y_l2)
            
            self.solution_l1 = yy_l1
            self.solution_l2 = yy_l2
            
            if self.iter_val==1 or self.iter_val == 2 or self.iter_val == 3 or self.iter_val == 4:
                x_l1 = "[ "
                x_l2 = "[ "
                for item in self.solution_l1:
                    x_l1 += f'{item:.{self.epsilon}f}' + " "
                x_l1 += "]"
                for item in self.solution_l2:
                    x_l2 += f'{item:.{self.epsilon}f}' + " "
                x_l2 += "]"

                print("In iteration ", self.iter_val, "\n\tfor lambda = 0.5 we have x = ", x_l1, "\n\tfor lambda = 0.9 we have x = ", x_l2)
                self.iter_val = self.iter_val + 1
            
            if norm(numpy.subtract(yy_l1, prev_sol_l1), ord = 2)< 0.00001 or norm(numpy.subtract(yy_l2, prev_sol_l2), ord = 2)< 0.00001:
                break
        x_l1 = "[ "
        x_l2 = "[ "
        for item in self.solution_l1:
            x_l1 += f'{item:.{self.epsilon}f}' + " "
        x_l1 += "]"
        for item in self.solution_l2:
            x_l2 += f'{item:.{self.epsilon}f}' + " "
        x_l2 += "]"
        print ("In the last iteration ", self.iter_val, "\n\tfor lambda = 0.5 we have x = " , x_l1, "\n\tfor lambda = 0.9 we have x = ", x_l2)
        res = 0
        for i in range(len(self.coefs)):
            res += self.coefs[i]*self.solution_l1[i]
        print("\nOptimum:", round(res, self.epsilon))

    def solve_max(self):
        self.iterior_point_method()

    def solve_min(self):
        for item in self.coefs:
            item *= -1
        self.iterior_point_method()

##
# function that will get data from user and check it for some typos, etc.
# and will envoke 2 methods (simplex, iterior-point) if exerything is fine =)
##
def userInput():
    type = input("Greetings, this programm will solve your LP problem using interior-point method.\nEnter the type of the problem(Max/Min): ").lower()
    if (type != "max" and type != "min"):
        print("ERROR: UNKNOWN TYPE")
        return
    print("NOTE: in input there shouldn't be any slack variables")
    try:
        objective_function = list(map(float, input("Enter the coefficients of the objective function: ").split(" ")))
        if(len(objective_function) == 0):
            print("ERROR: NO COEFFICIENTS")
            return
        
        amount = int(input("Enter amount of the constraints(not assuming x>=0): "))
        if(amount < 1):
            print("ERROR: AMOUNT < 1 ?!")
        constraints = []
        for i in range(amount):
            constraint = list(map(float, input(f"Enter the {i+1} constraint function coefficients: ").split(" ")))
            if(len(constraint) == 0):
                print("ERROR: NOT ENOUGH COEFFICIENTS")
                return
            constraints.append(constraint)
            
        right_hand_side = list(map(float, input("Enter the right-hand side numbers: ").split(" ")))
        if(len(right_hand_side) != amount):
            print("ERROR: NOT ENOUGH COEFFICIENTS")
            return
        else:
            for i in right_hand_side:
                if(i < 0):
                    print("The method is not applicable!")
                    return

        accuracy = input("Enter the approximation accuracy: ")

        if int(accuracy) > 0: 
            accuracy = int(accuracy)
        else:
            k = 0
            accuracy = float(accuracy)
            while True:
                k += 1
                accuracy *= 10
                if accuracy > 0:
                    break
            accuracy = k

        lp_simplex = Simplex(objective_function, constraints, right_hand_side, accuracy)

        #the following block of code make change of input data for different methods
        for i in range(len(constraints)):
            objective_function.append(0)
            constraint_iter = constraints[i]
            for j in range(len(constraints)):
                if i == j:
                    constraint_iter.append(1)
                else:
                    constraint_iter.append(0)

        lp_iterior_point = IteriorPoint(objective_function, constraints, right_hand_side, accuracy)
        # lp_simplex.print_initial()

        if(type == "max"):
            print("===============================================================")
            print ("SIMPLEX")
            lp_simplex.solve_maximize()
            print("===============================================================")
            print ("INTERIOR POINT")
            lp_iterior_point.solve_max()
            print("===============================================================")


        else:
            print("===============================================================")
            print ("SIMPLEX")
            lp_simplex.solve_minimize()
            print("===============================================================")
            print ("INTERIOR POINT")
            lp_iterior_point.solve_min()
            print("===============================================================")
        

    except ValueError:
        print("ERROR: NOT A NUMBER")
        return

userInput()

# initial_solution = numpy.array([1, 1, 1, 315, 174, 169], float)
# constraints = numpy.array([[18, 15, 12, 1, 0, 0], [6, 4, 8, 0, 1, 0], [5, 3, 3, 0, 0, 1]], float)
# coefs = numpy.array([9, 10, 16, 0, 0, 0], float)

# lp_simplex = Simplex([9,10,16], [[18, 15, 12], [6, 4, 8], [5, 3, 3]], [360, 192, 180], 4)
# lp_iterior_point = IteriorPoint(coefs, constraints, [360, 192, 180], 10)
# lp_iterior_point.solve_max()

# lp_simplex.solve_maximize()
# # lp_simplex.solve_minimize()
# lp_simplex.print_solved()
