function [] = rmsdelayspread(Times,CBB)
                
            %Calculate power-delay profile and various metrics
            df  = (Frequencies(2) - Frequencies(1)) * 1000000000.0; dt = (Times(2, 1) - Times(1, 1)) * 0.000000001;
            PNormPDPscale = sqrt(Filter2Avg / (N * df)) / FilterAvg;

%            Dim DelayAvg As Double = 0.0, PowerThreshTotal As Double = 0.0
            tau_a = 0.0; tau_a_fit = 0.0; tau_a_index = 0; FirstPeakValue = 0.0; 
            tau_a_max = 0.0; tau_a_max_fit = 0.0; tau_a_max_index = 0;   %tau_a defined as maximum peak
            Threshhold = 10.^(-abs(Threshold_dB) / 10.0);
            DIThreshhold = 10.^(-abs(DITThreshold_dB) / 10.0);
            MCThreshhold = 10.0.^(-abs(MCThreshhold_dB)/10.0);

            %Find the peak power
            ndx1 = find(abs(CBB) == max(abs(CBB)));
            MCmax = abs(CBB(ndx1)).^2;
            tau_a_max = times(ndx1);
            %for kT As Integer = 1 To Times.NRows
            %    Dim Amplitude As Double = Abs(CBB(kT))  %Get both the power and thresholded power
            %    Dim Power As Double = Amplitude * Amplitude
            %    If MCmax < Power Then MCmax = Power : tau_a_max = Times(kT) : tau_a_max_index = kT 'The max power
%             Next kT

            %Define the Multipath Component threshold with respect to the peak component.
            MCThreshhold = MCmax.* MCThreshhold;                         %Define multipath threshold with respect to peak power
            if (MCThreshhold < Threshhold); MCThreshhold = Threshhold; end  %Set to noise threshold if it is below the noise threshold

            %Step through the times and find the average and thresholded powers, and the avearage delay.
            timelength = length(times);
            for kT = 1:timelength   %Step through the times
                Amplitude = abs(CBB(kT));  %Get both the power and thresholded power
                Power = Amplitude.*Amplitude;
                PowerTresh = Power;     %The thresh-holded power
                if (Power < Threshhold); PowerTresh = 0.0; end
                PowerThreshed(kT) = PowerTresh;      %Keep track of the thresh-holded power for later
                %'Get the start of the delay interval, look for first and last change in sign
                if kT > 1
                    if (PowerLast - DIThreshhold) * (Power - DIThreshhold) <= 0.0; %Found a change in sign, get intercept
                        DI2 = (Times(kT - 1) + Times(kT)) / 2.0;
                        if Power ~= PowerLast; DI2 = Times(kT - 1) + (DIThreshhold - PowerLast) * (Times(kT) - Times(kT - 1)) / (Power - PowerLast);end
                        if Not DIfound;
                            DI1 = DI2;   %The first chang in sign
                        end
                        DIfound = True;
                        end
                    end
                %Generate the amplitudes used to create the power-delay profile
                PDP(kT, 1) = Times(kT) ; PDP(kT, 2) = Amplitude;
                %Normalized PDP so that integrating PDP of a channel with S21=1 wrt time gives 1. Implies that multiplying by power gives actual PDP, with int PDP dt = power in time interval.
                PNormPDP(kT, 1) = Times(kT) ; PNormPDP(kT, 2) = Amplitude / PNormPDPscale;
                %Normalized PDP so that summing elements of PDP of a channel with S21=1 gives 1.
                DiscretePNormPDP(kT, 1) = Times(kT) ; DiscretePNormPDP(kT, 2) = PNormPDP(kT, 2) * Math.Sqrt((Times(2) - Times(1)) * 0.000000001);
                %Rescale the normalized PDP by 1/sqrt(BWareascale)
                PNormPDP(kT, 2) = PNormPDP(kT, 2) / sqrt(BWareascale);
                %Generate the cumulitive power-delay profile
                if kT == 1;
                    PDP_CD(kT) = PowerTresh;
                else
                    PDP_CD(kT) = PDP_CD(kT - 1) + PowerTresh;
                end

                %Sum up quantities needed to calculate average delay
                DelayAvg += PowerTresh * Times(kT); PowerThreshTotal += PowerTresh;

                %tau_a defined in ITU-% P.1407-5, page 5, section 2.2.2
                %Dim tau_a As Double = 0.0, FirstPeakValue As Double = 0.0, FirstPeak As Boolean = True
                if FirstPeak                            %Make sure we are still on or before the first peak
                    if Power >= MCThreshhold    %Make sure we are on the first peak
                        if FirstPeakValue <= Power      %We are climbing the first peak
                            FirstPeakValue = Power ;         %Keep track of where we are
                            tau_a = Times(kT);
                            tau_a_index = kT;
                        else                                %We just started to descend off of the first peak
                            FirstPeak = False;               %Stop looking for the first peak
                        end
                    else
                        if FirstPeakValue > 0.0; FirstPeak = False; end  %We were on the peak, but the power dropped to zero. So we must have descended from the peak
                    end
                end

                %Keep track of the power we had last time in the calculation
                PowerLast = Power;
            end
            DelayAvg = DelayAvg./PowerThreshTotal;  %The power-weighted delay in absolute time

            %find the number of multipath componenets
            MCnumber = 0.0;
            for kT = 1:timelength-1
                if (PowerThreshed(kT) - MCThreshhold) * (PowerThreshed(kT + 1) - MCThreshhold) <= 0.0; MCnumber += 0.5; end
            end

            %Find the delay window
            Dim DWThreshhold As Double = MechValues(3).MechanismValue(MechanismList1)   
            If DWThreshhold < 0.5 Or DWThreshhold >= 1.0 Then DWThreshhold = 0.9
            Dim DWT1 As Double = (1.0 - DWThreshhold) / 2.0, DWT2 As Double = 1.0 - DWT1
            Dim kT_5 As Integer = 0, kT_95 As Integer = 0
            for kT = 1:timelength-1    %Step through the times
                If (PDP_CD(kT) - DWT1 * PDP_CD(Times.NRows)) * (PDP_CD(kT + 1) - DWT1 * PDP_CD(Times.NRows)) <= 0 Then kT_5 = kT
                If (PDP_CD(kT) - DWT2 * PDP_CD(Times.NRows)) * (PDP_CD(kT + 1) - DWT2 * PDP_CD(Times.NRows)) <= 0 Then kT_95 = kT
            end
            Dim t_5 As Double = 0.0, t_95 As Double = 0.0   'Default valuse if we could not find any threshhold crossings
            If kT_5 > 0 And kT_95 > 0 Then                   'We found some threshhold crossings
                If kT_5 = 1 Then kT_5 = 2
                t_5 = (Times(kT_5 + 1) + Times(kT_5)) / 2.0
                If PDP_CD(kT_5 + 1) <> PDP_CD(kT_5) Then t_5 = Times(kT_5) + (DWT1 * PDP_CD(Times.NRows) - PDP_CD(kT_5)) * (Times(kT_5 + 1) - Times(kT_5)) / (PDP_CD(kT_5 + 1) - PDP_CD(kT_5))
                t_95 = (Times(kT_95 + 1) + Times(kT_95)) / 2.0
                If PDP_CD(kT_95 + 1) <> PDP_CD(kT_95) Then t_95 = Times(kT_95) + (DWT2 * PDP_CD(Times.NRows) - PDP_CD(kT_95)) * (Times(kT_95 + 1) - Times(kT_95)) / (PDP_CD(kT_95 + 1) - PDP_CD(kT_95))
            End If

            'Check that PNormPDP is calculated correctly for S21=1
            Dim Debug As Boolean = False
            If Debug Then
                Dim Sum As Double = 0.0, Sum1 As Double = 0.0, Sum2 As Double = 0.0, Second As Boolean = True
                For kT As Integer = 1 To Times.NRows    'Step through the times
                    Sum += (PNormPDP(kT, 2) ^ 2) * dt
                    Sum1 += (PDP(kT, 2) ^ 2) * dt
                    If Second Then Sum2 += DiscretePNormPDP(kT, 2) ^ 2
                    'Second = Not Second
                Next kT
            End If

            'Find rms delay spread
            Dim RMSDelaySpread As Double = 0.0
            For kT As Integer = 1 To Times.NRows    'Step through the times
                RMSDelaySpread += PowerThreshed(kT) * ((Times(kT) - DelayAvg) ^ 2)
            Next kT
            RMSDelaySpread = Math.Sqrt(RMSDelaySpread / PowerThreshTotal)

return
% 
%         ''' <summary>
%         ''' Fit the arrival time of the first peak
%         ''' </summary>
%         ''' <param name="PeakIndex">Peak arrival time index from ITU definition</param>
%         ''' <param name="Times">An array of the times.</param>
%         ''' <param name="PowerThreshed">An array of the threshholded powers.</param>
%         ''' <returns>This is much better than choosing the maximum value when picking the arrival time, as done in the ITU.
%         ''' Choosing the maximum value implies integer time points, and uncertainties may not make much sense.</returns>
%         Private Function PeakCenter(ByVal PeakIndex As Integer, ByRef Times As RealMatrix, ByRef PowerThreshed As RealMatrix) As Double
% 
%             'Get a reasonable array of times around the peak.
%             Dim Kfirst As Integer = PeakIndex, Klast As Integer = PeakIndex, PeakMax As Double = PowerThreshed(PeakIndex)
%             While Kfirst > 1 And PowerThreshed(Kfirst) / PeakMax > 0.3
%                 Kfirst += -1
%             End While
%             While Klast < PowerThreshed.NRows And PowerThreshed(Klast) / PeakMax > 0.3
%                 Klast += 1
%             End While
%             Dim PowerThreshedNearPeak As New RealMatrix(1 + Klast - Kfirst), TimesNearPeak As New RealMatrix(1 + Klast - Kfirst)
%             For kT As Integer = 1 To PowerThreshedNearPeak.NRows
%                 PowerThreshedNearPeak(kT) = PowerThreshed((Kfirst - 1) + kT)
%                 TimesNearPeak(kT) = Times((Kfirst - 1) + kT)
%             Next kT
% 
%             'Pick a reasonable number of points to use to get the rising and falling edges
%             Dim nfit As Integer = PowerThreshedNearPeak.NRows / 4
%             Dim myOrder As Integer = 0
%             If nfit < 3 Then myOrder = 1
% 
%             'Find the maximimum
%             If nfit < 2 Then    'Not enough data to fit the rising and falling edges. Punt.
%                 Return Times(PeakIndex)
%             Else
%                 Dim Tr As Double = findtrImp_wls2(TimesNearPeak, PowerThreshedNearPeak, 0.5, nfit, -1, myOrder, 0, False)  'Rising edge
%                 Dim Tf As Double = findtrImp_wls2(TimesNearPeak, PowerThreshedNearPeak, 0.5, nfit, 1, myOrder, 0, False)   'Falling edge
%                 Return 0.5 * (Tr + Tf)  'The center of the peak
%             End If
% 
%         End Function
% 
%     End Class
