Reference: ''Shuo Li, Lixia Xiao, and Tao Jiang. An Efficient Matching Pursuit Based Compressive Sensing Detector For Uplink Grant-Free NOMA. in IEEE Transactions on Vehicular Technology, vol. 70, no. 2, pp:2012-2017, Feb. 2021.''

The implementation of block CoSaMP and block OMP is partially based on an open-source project for CoSaMP and OMP.
https://github.com/rasikraj01/CompressiveSensing/tree/399e35b8921d08d787c23c52576ab42331ebd66d

Note that the LS estimation is not used in block CoSaMP, since the CS_CoSaMP does not  applied it. However, LS estimation can be added for BCoSaMP to improve the BER performance.